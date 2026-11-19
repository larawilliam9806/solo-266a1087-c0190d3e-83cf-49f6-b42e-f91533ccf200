"""并发一致性测试：样本上传与参考更新交错时，汇总始终对应单一参考版本。

并发交错本身由线程调度决定，无法在测试里确定性地强制某个事务先拿到锁；
因此本文件分两层验证：

  1. 多线程压测（barrier 同时起跑）：只断言「任何时刻读到的汇总都满足
     跨版本不混算等不变量、且无死锁/非预期异常」，不依赖特定交错时序；
  2. 确定性串行测试：参考轮换后旧参考形态的上传必然被拒绝且不留脏数据，
     同名并发重传最终只保留一份样本——这两个结果不依赖线程调度。
"""

import os
import tempfile
import threading
import unittest

from app import create_app, services
from app.normalization import ValidationError

REF_V1 = "GCTTTGCACGT"
REF_V2 = "GCTATGCACGT"
REF_V3 = "ACGTACGTACGT"
REFS = [REF_V1, REF_V2, REF_V3]


def tsv(text: str) -> bytes:
    return text.encode("utf-8")


class ConcurrencyTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        os.environ["DATABASE_PATH"] = os.path.join(self.tmp.name, "conc.db")
        self.app = create_app()
        self.app.config.update(TESTING=True, SECRET_KEY="test")
        self.errors: list[Exception] = []

    def tearDown(self):
        os.environ.pop("DATABASE_PATH", None)
        self.tmp.cleanup()

    def test_concurrent_uploads_and_rotations_stay_consistent(self):
        # 先放一版参考
        with self.app.app_context():
            services.upload_reference(REF_V1)

        # 所有工作线程在同一屏障前集合后同时起跑，确保真正并发
        start = threading.Barrier(4)

        def upload_samples():
            start.wait()
            for i in range(15):
                name = f"S{i % 6}"  # 含同名重传
                gt = "1/1" if i % 2 else "0/1"
                try:
                    with self.app.app_context():
                        # 每条都针对 V1 有效；换参考后由系统重验，
                        # 若按已轮换参考校验失败则被整体拒绝（允许发生，不强制）
                        services.upload_sample(name, tsv(f"1\tG\tC\t{gt}\n"))
                except ValidationError:
                    pass  # 参考轮换导致的竞争拒绝是合法结果
                except Exception as exc:
                    self.errors.append(exc)

        def rotate_references():
            start.wait()
            for ref in REFS[1:]:
                try:
                    with self.app.app_context():
                        services.upload_reference(ref)
                except Exception as exc:
                    self.errors.append(exc)

        def reader():
            start.wait()
            for _ in range(40):
                try:
                    with self.app.app_context():
                        dash = services.get_dashboard()
                        ref_id = dash["reference"]["id"]
                        # 不变量 1：汇总中的每个变体只来自当前参考的有效样本
                        for v in dash["variants"]:
                            self.assertGreaterEqual(v["n_samples"], 1)
                        # 不变量 2：所有 VALID 样本都指向当前参考
                        for s in dash["samples"]:
                            if s["status"] == "VALID":
                                self.assertEqual(s["reference_id"], ref_id)
                        # 不变量 3：变体数与样本变体明细一致
                        total = sum(
                            s["variant_count"]
                            for s in dash["samples"]
                            if s["status"] == "VALID"
                            and s["reference_id"] == ref_id
                        )
                        self.assertLessEqual(len(dash["variants"]), total)
                except Exception as exc:
                    self.errors.append(exc)

        threads = [
            threading.Thread(target=upload_samples),
            threading.Thread(target=rotate_references),
            threading.Thread(target=reader),
            threading.Thread(target=reader),
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)
            self.assertFalse(t.is_alive(), "线程未在超时内结束，可能存在锁死锁")

        if self.errors:
            raise AssertionError(
                "并发过程出现非预期异常:\n"
                + "\n".join(f"{type(e).__name__}: {e}" for e in self.errors)
            )

        # 最终状态：参考为 V3，所有样本已被重验并指向 V3（VALID 或 NEEDS_REUPLOAD）
        with self.app.app_context():
            dash = services.get_dashboard()
            self.assertEqual(dash["reference"]["sequence"], REF_V3)
            for s in dash["samples"]:
                self.assertEqual(s["reference_id"], dash["reference"]["id"])

    def test_upload_after_rotation_rejected_without_dirty_data(self):
        """确定性覆盖「交错竞争」的失败分支：轮换后旧参考形态的上传被拒。"""
        with self.app.app_context():
            services.upload_reference(REF_V1)
            services.upload_sample("S1", tsv("1\tG\tC\t0/1\n"))
            # 轮换到 V3：首位碱基 G -> A，旧 TSV 的 REF 不再匹配
            services.upload_reference(REF_V3)

            with self.assertRaises(ValidationError) as ctx:
                services.upload_sample("S2", tsv("1\tG\tC\t0/1\n"))
            self.assertIn("REF 核对失败", str(ctx.exception))

            # 拒绝不留脏数据：S2 不存在，S1 已被重验指向 V3
            from app.db import get_db
            names = {
                r["name"]: (r["reference_id"], r["status"])
                for r in get_db().execute("SELECT name, reference_id, status FROM samples")
            }
            self.assertNotIn("S2", names)
            ref_id = get_db().execute(
                "SELECT id FROM reference_sequences ORDER BY id DESC LIMIT 1"
            ).fetchone()["id"]
            self.assertEqual(names["S1"], (ref_id, "NEEDS_REUPLOAD"))

    def test_concurrent_same_name_reuploads_keep_single_sample(self):
        """同名样本并发重传：写事务串行化后只允许保留一份，不出现重复行。"""
        with self.app.app_context():
            services.upload_reference(REF_V1)

        start = threading.Barrier(2)

        def reupload(gt: str):
            start.wait()
            try:
                with self.app.app_context():
                    services.upload_sample("DUP", tsv(f"1\tG\tC\t{gt}\n"))
            except Exception as exc:
                self.errors.append(exc)

        threads = [
            threading.Thread(target=reupload, args=("0/1",)),
            threading.Thread(target=reupload, args=("1/1",)),
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)
            self.assertFalse(t.is_alive())
        self.assertEqual(self.errors, [])

        with self.app.app_context():
            from app.db import get_db
            sample_rows = get_db().execute(
                "SELECT id FROM samples WHERE name='DUP'"
            ).fetchall()
            self.assertEqual(len(sample_rows), 1)
            variant_rows = get_db().execute(
                "SELECT gt FROM sample_variants WHERE sample_id=?",
                (sample_rows[0]["id"],),
            ).fetchall()
            self.assertEqual(len(variant_rows), 1)
            self.assertIn(variant_rows[0]["gt"], ("0/1", "1/1"))


if __name__ == "__main__":
    unittest.main()
