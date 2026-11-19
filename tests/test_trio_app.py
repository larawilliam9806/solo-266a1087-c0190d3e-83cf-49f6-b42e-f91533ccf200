"""三联核查端到端测试：Flask 测试客户端 + 临时 SQLite。"""

import io
import os
import tempfile
import unittest

from app import create_app, services
from app.normalization import ValidationError

# 12 bp 参考；SNP 位点 1/2/3/5/8/10/11，indel 位点 4
REF_V1 = "GCTTTGCACGTA"
# 第 4 位 T->A，全部 SNP 的 REF 不变，indel (3,T,TT) 的锚碱基也不变，
# 用于「参考重验后再次核查」场景
REF_V2 = "GCTATGCACGTA"


def tsv(text: str) -> bytes:
    return text.encode("utf-8")


class TrioAppTestBase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        os.environ["DATABASE_PATH"] = os.path.join(self.tmp.name, "trio.db")
        self.app = create_app()
        self.app.config.update(TESTING=True, SECRET_KEY="test")
        self.client = self.app.test_client()
        with self.app.app_context():
            services.upload_reference(REF_V1)

    def tearDown(self):
        os.environ.pop("DATABASE_PATH", None)
        self.tmp.cleanup()

    def upload(self, name: str, content: str):
        data = {"name": name, "sample": (io.BytesIO(tsv(content)), f"{name}.tsv")}
        return self.client.post(
            "/samples", data=data, content_type="multipart/form-data",
            follow_redirects=True,
        )

    def sample_id(self, name: str) -> int:
        with self.app.app_context():
            from app.db import get_db
            return get_db().execute(
                "SELECT id FROM samples WHERE name=?", (name,)
            ).fetchone()["id"]

    def check(self, father: str, mother: str, child: str):
        resp = self.client.post(
            "/trio-check",
            data={
                "father": str(self.sample_id(father)),
                "mother": str(self.sample_id(mother)),
                "child": str(self.sample_id(child)),
            },
            follow_redirects=True,
        )
        return resp

    def check_result(self, father: str, mother: str, child: str) -> dict:
        with self.app.app_context():
            return services.run_trio_check(
                {
                    "father": self.sample_id(father),
                    "mother": self.sample_id(mother),
                    "child": self.sample_id(child),
                }
            )


class TrioCheckWorkflowTest(TrioAppTestBase):
    def setUp(self):
        super().setUp()
        # 父亲：pos1 检出 C(0/1)、pos2 检出 T(1/1)、pos3 0/0（已检测未检出 T）、
        #       pos5 检出 G(1/1)、pos8 检出 C(0/1)、外加一条插入（核查应忽略 indel）
        self.upload(
            "Dad",
            "POS\tREF\tALT\tGT\n"
            "1\tG\tC\t0/1\n2\tC\tT\t1/1\n3\tT\tA\t0/0\n"
            "5\tT\tG\t1/1\n8\tA\tC\t0/1\n3\tT\tTT\t0/1\n",
        )
        # 母亲：pos1 0/0、pos2 0/0、pos3 0/0、pos5 0/0、pos10 检出 T(0/1)
        self.upload(
            "Mom",
            "1\tG\tC\t0/0\n2\tC\tT\t0/0\n3\tT\tA\t0/0\n"
            "5\tT\tG\t0/0\n10\tG\tT\t0/1\n",
        )
        # 子代：
        #   pos1 0/0  -> 父 0/1 母 0/0 可解释
        #   pos2 1/1  -> 父 1/1 母 0/0 不可解释（母亲传不出 ALT）-> 违规
        #   pos3 0/0  -> 三方 0/0 可解释（纯 0/0 位点也参与核查）
        #   pos5 0/1  -> 父 1/1 母 0/0 可解释
        #   pos8 0/1  -> 母亲缺记录 -> 排除
        #   pos10 1/1 -> 父亲缺记录 -> 排除
        self.upload(
            "Kid",
            "1\tG\tC\t0/0\n2\tC\tT\t1/1\n3\tT\tA\t0/0\n"
            "5\tT\tG\t0/1\n8\tA\tC\t0/1\n10\tG\tT\t1/1\n",
        )

    def test_full_classification(self):
        result = self.check_result("Dad", "Mom", "Kid")
        self.assertEqual(
            sorted(v["pos"] for v in result["consistent"]), [1, 3, 5]
        )
        self.assertEqual([v["pos"] for v in result["violations"]], [2])
        self.assertEqual(
            sorted(v["pos"] for v in result["excluded"]), [8, 10]
        )
        # 总位点数只统计 SNP（6 个），indel 不计入
        self.assertEqual(result["n_loci"], 6)
        viol = result["violations"][0]
        self.assertEqual(viol["alts"], ["T"])
        gts = {r["role"]: r["gt"] for r in viol["rows"]}
        self.assertEqual(gts, {"father": "1/1", "mother": "0/0", "child": "1/1"})
        for item in result["excluded"]:
            self.assertIn("三方记录不齐全", item["reason"])

    def test_result_page_renders_three_sections(self):
        resp = self.check("Dad", "Mom", "Kid")
        self.assertEqual(resp.status_code, 200)
        html = resp.get_data(as_text=True)
        self.assertIn("无法由双亲各传一个等位基因解释", html)
        self.assertIn("缺记录或出现多个 ALT", html)
        self.assertIn("遗传可解释的位点", html)
        # 违规位点 pos2 与三方 GT 均展示
        self.assertIn(">2<", html)
        self.assertIn("Dad", html)
        self.assertIn("Mom", html)
        self.assertIn("Kid", html)
        # 参考版本与核查时刻展示
        self.assertIn("参考 #1", html)
        self.assertIn("核查时刻", html)

    def test_zero_zero_only_sample_participates(self):
        # 新增一份只含 0/0 记录的样本，仍可被选为核查方
        self.upload(
            "Zero",
            "1\tG\tC\t0/0\n2\tC\tT\t0/0\n3\tT\tA\t0/0\n",
        )
        result = self.check_result("Dad", "Mom", "Zero")
        pos1 = next(v for v in result["consistent"] if v["pos"] == 1)
        child_row = next(r for r in pos1["rows"] if r["role"] == "child")
        self.assertEqual(child_row["gt"], "0/0")
        # 该样本在首页汇总里没有任何检出贡献
        page = self.client.get("/").get_data(as_text=True)
        self.assertIn("Zero", page)

    def test_dashboard_and_detail_exclude_zero_zero(self):
        with self.app.app_context():
            dash = services.get_dashboard()
        # pos3 的三方记录都是 0/0（Dad 另有同 ALT 的 0/0），不应出现在汇总
        keys = {(v["pos"], v["ref"], v["alt"]) for v in dash["variants"]}
        self.assertNotIn((3, "T", "A"), keys)
        self.assertIn((1, "G", "C"), keys)  # Dad 0/1 检出
        # 详情页只有检出者 Dad，Mom/Kid 的 0/0 不展示
        detail = self.client.get("/variants/1/G>C").get_data(as_text=True)
        self.assertIn("Dad", detail)
        self.assertNotIn("Mom", detail)
        self.assertNotIn("Kid", detail)

    def test_multi_alt_locus_excluded(self):
        # 让母亲在 pos1 申报另一种 ALT（多等位基因位点）
        self.upload(
            "Mom",
            "1\tG\tA\t0/1\n2\tC\tT\t0/0\n3\tT\tA\t0/0\n"
            "5\tT\tG\t0/0\n10\tG\tT\t0/1\n",
        )
        result = self.check_result("Dad", "Mom", "Kid")
        pos1 = next(v for v in result["excluded"] if v["pos"] == 1)
        self.assertEqual(pos1["alts"], ["A", "C"])
        self.assertIn("多等位基因", pos1["reason"])
        self.assertNotIn(1, [v["pos"] for v in result["violations"]])

    def test_duplicate_selection_rejected(self):
        sid = self.sample_id("Dad")
        with self.app.app_context():
            with self.assertRaises(ValidationError) as ctx:
                services.run_trio_check(
                    {"father": sid, "mother": sid, "child": self.sample_id("Kid")}
                )
        self.assertIn("三个不同的样本", str(ctx.exception))
        # 路由层：失败后回到首页并提示重选
        resp = self.client.post(
            "/trio-check",
            data={"father": sid, "mother": sid, "child": sid},
            follow_redirects=True,
        )
        self.assertIn("请重新选择", resp.get_data(as_text=True))

    def test_stale_sample_id_prompts_reselect(self):
        sid = self.sample_id("Mom")
        # 删除样本后旧 id 失效
        self.client.post(f"/samples/{sid}/delete", follow_redirects=True)
        with self.app.app_context():
            with self.assertRaises(ValidationError) as ctx:
                services.run_trio_check(
                    {
                        "father": self.sample_id("Dad"),
                        "mother": sid,
                        "child": self.sample_id("Kid"),
                    }
                )
        self.assertIn("已失效", str(ctx.exception))

    def test_revalidation_invalidates_then_reselect_uses_current_data(self):
        # 第一次核查可行
        before = self.check_result("Dad", "Mom", "Kid")
        self.assertEqual([v["pos"] for v in before["violations"]], [2])

        # 参考重验：V2 下三份样本的 SNP 与 indel 锚碱基仍匹配，全部保持 VALID，
        # 再次核查使用当前参考的数据
        with self.app.app_context():
            services.upload_reference(REF_V2)
        after = self.check_result("Dad", "Mom", "Kid")
        self.assertEqual(after["reference"]["id"], before["reference"]["id"] + 1)
        self.assertEqual(
            sorted(v["pos"] for v in after["consistent"]), [1, 3, 5]
        )

        # 令 Mom 重传为与 V2 不匹配的文件 -> 整份上传失败，原核查依据保留
        bad = self.upload("Mom", "1\tT\tC\t0/0\n")  # pos1 参考实为 G
        self.assertIn("REF 核对失败", bad.get_data(as_text=True))
        unchanged = self.check_result("Dad", "Mom", "Kid")
        self.assertEqual(
            sorted(v["pos"] for v in unchanged["consistent"]), [1, 3, 5]
        )

        # Mom 重传为有效新内容（去掉 pos2 记录）后，再次核查自动使用当前数据：
        # pos2 变为「母亲缺记录」进入排除清单，违规清单清空
        self.upload(
            "Mom",
            "1\tG\tC\t0/0\n3\tT\tA\t0/0\n5\tT\tG\t0/0\n10\tG\tT\t0/1\n",
        )
        renewed = self.check_result("Dad", "Mom", "Kid")
        self.assertEqual(renewed["violations"], [])
        pos2 = next(v for v in renewed["excluded"] if v["pos"] == 2)
        self.assertIn("母亲「Mom」", pos2["reason"])

    def test_needs_reupload_sample_rejected(self):
        # 上传一份在 V2 下会 REF 不符的样本，换参考后变待重传
        self.upload("Fragile", "4\tT\tTT\t1/1\n")
        with self.app.app_context():
            services.upload_reference(REF_V2)
        with self.app.app_context():
            from app.db import get_db
            row = get_db().execute(
                "SELECT id, status FROM samples WHERE name='Fragile'"
            ).fetchone()
            self.assertEqual(row["status"], "NEEDS_REUPLOAD")
            with self.assertRaises(ValidationError) as ctx:
                services.run_trio_check(
                    {
                        "father": row["id"],
                        "mother": self.sample_id("Mom"),
                        "child": self.sample_id("Kid"),
                    }
                )
        self.assertIn("不是当前参考", str(ctx.exception))


class TrioFormAvailabilityTest(TrioAppTestBase):
    def test_form_hidden_until_three_valid_samples(self):
        page = self.client.get("/").get_data(as_text=True)
        self.assertIn("有效样本不足 3 个", page)
        self.upload("S1", "1\tG\tC\t0/1\n")
        self.upload("S2", "1\tG\tC\t1/1\n")
        page = self.client.get("/").get_data(as_text=True)
        self.assertIn("有效样本不足 3 个", page)
        self.upload("S3", "1\tG\tC\t0/0\n")  # 仅 0/0 也是有效样本
        page = self.client.get("/").get_data(as_text=True)
        self.assertIn("三联核查", page)
        self.assertIn('name="father"', page)
        self.assertIn('name="mother"', page)
        self.assertIn('name="child"', page)
        self.assertIn("S3", page)

    def test_check_without_reference_rejected(self):
        # 清空参考场景：直接调用服务应给出明确错误
        with self.app.app_context():
            from app.db import get_db
            get_db().execute("DELETE FROM reference_sequences")
            with self.assertRaises(ValidationError) as ctx:
                services.run_trio_check({"father": 1, "mother": 2, "child": 3})
        self.assertIn("尚未上传参考", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
