"""业务逻辑：参考/样本上传、重验、替换、汇总查询与三联核查。"""

from __future__ import annotations

import sqlite3
from datetime import datetime, timezone

from .db import read_tx, write_tx
from .normalization import (
    MAX_NAME_LEN,
    ValidationError,
    parse_reference,
    parse_tsv,
    validate_and_normalize,
)
from .trio import ROLES, summarize_trio_loci


# --------------------------------------------------------------------------- #
# 写操作
# --------------------------------------------------------------------------- #
def upload_reference(raw: bytes) -> dict:
    """上传新参考；若已存在样本，全部按新参考重新校验。

    成功后：
      * 通过重验的样本 -> VALID，变体替换为新参考下的归一化结果；
      * REF 不符 / 归一化后重复的样本 -> NEEDS_REUPLOAD，变体删除，记录原因，
        从汇总剔除；
      * 旧参考及旧变体在同一事务中整体替换，失败则全部回滚。
    """
    sequence = parse_reference(raw)
    with write_tx() as conn:
        cur = conn.execute("SELECT id FROM samples")
        old_samples = cur.fetchall()

        cur = conn.execute(
            "INSERT INTO reference_sequences (sequence, length) VALUES (?, ?)",
            (sequence, len(sequence)),
        )
        new_ref_id = cur.lastrowid

        revalidated, failed = [], []
        for row in old_samples:
            outcome = _revalidate_sample(conn, row["id"], new_ref_id, sequence)
            (revalidated if outcome else failed).append(
                conn.execute(
                    "SELECT name, status, error FROM samples WHERE id=?",
                    (row["id"],),
                ).fetchone()
            )
    return {
        "reference_id": new_ref_id,
        "length": len(sequence),
        "revalidated": revalidated,
        "failed": failed,
    }


def upload_sample(name: str, raw: bytes) -> dict:
    """上传（或同名重传）一份样本 TSV。

    重传时先完整校验新文件，成功才替换旧版本；旧变体在同一事务中删除、
    新变体写入，任何错误整体回滚，当前结果不受影响。
    """
    name = (name or "").strip()
    if not name:
        raise ValidationError("样本名称不能为空。")
    if len(name) > MAX_NAME_LEN:
        raise ValidationError(f"样本名称不能超过 {MAX_NAME_LEN} 个字符。")

    # 解析在事务外完成（不触碰数据库，不影响现有结果）
    rows = parse_tsv(raw)
    raw_text = raw.decode("utf-8")

    with write_tx() as conn:
        ref_row = conn.execute(
            "SELECT id, sequence FROM reference_sequences ORDER BY id DESC LIMIT 1"
        ).fetchone()
        if ref_row is None:
            raise ValidationError("请先上传参考序列，再上传样本。")
        ref_id, ref_seq = ref_row["id"], ref_row["sequence"]

        # 校验 + 左对齐 + 去冗余；失败抛 ValidationError，事务回滚
        variants = validate_and_normalize(rows, ref_seq)

        old = conn.execute("SELECT id FROM samples WHERE name=?", (name,)).fetchone()
        if old is not None:
            conn.execute("DELETE FROM sample_variants WHERE sample_id=?", (old["id"],))
            conn.execute("DELETE FROM samples WHERE id=?", (old["id"],))

        cur = conn.execute(
            "INSERT INTO samples (name, reference_id, status, error, raw_tsv) "
            "VALUES (?, ?, 'VALID', NULL, ?)",
            (name, ref_id, raw_text),
        )
        sample_id = cur.lastrowid
        _insert_variants(conn, sample_id, variants)

    return {"sample": name, "count": len(variants), "replaced": old is not None}


def delete_sample(sample_id: int) -> None:
    with write_tx() as conn:
        conn.execute("DELETE FROM samples WHERE id=?", (sample_id,))


def _revalidate_sample(
    conn: sqlite3.Connection, sample_id: int, ref_id: int, ref_seq: str
) -> bool:
    """按给定参考重新校验一个样本，返回是否通过。"""
    row = conn.execute(
        "SELECT raw_tsv FROM samples WHERE id=?", (sample_id,)
    ).fetchone()
    try:
        rows = parse_tsv(row["raw_tsv"])
        variants = validate_and_normalize(rows, ref_seq)
    except ValidationError as exc:
        conn.execute("DELETE FROM sample_variants WHERE sample_id=?", (sample_id,))
        conn.execute(
            "UPDATE samples SET reference_id=?, status='NEEDS_REUPLOAD', "
            "error=?, revalidated_at=datetime('now') WHERE id=?",
            (ref_id, str(exc), sample_id),
        )
        return False

    conn.execute("DELETE FROM sample_variants WHERE sample_id=?", (sample_id,))
    conn.execute(
        "UPDATE samples SET reference_id=?, status='VALID', error=NULL, "
        "revalidated_at=datetime('now') WHERE id=?",
        (ref_id, sample_id),
    )
    _insert_variants(conn, sample_id, variants)
    return True


def _insert_variants(
    conn: sqlite3.Connection, sample_id: int, variants: list[dict]
) -> None:
    conn.executemany(
        "INSERT INTO sample_variants "
        "(sample_id, pos, ref, alt, gt, kind, orig_pos, orig_ref, orig_alt) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            (
                sample_id,
                v["pos"],
                v["ref"],
                v["alt"],
                v["gt"],
                v["kind"],
                v["orig_pos"],
                v["orig_ref"],
                v["orig_alt"],
            )
            for v in variants
        ],
    )


# --------------------------------------------------------------------------- #
# 读操作（同一读事务/快照内完成，页面只展示同一版参考的完整汇总）
# --------------------------------------------------------------------------- #
def get_dashboard() -> dict:
    with read_tx() as conn:
        ref = conn.execute(
            "SELECT * FROM reference_sequences ORDER BY id DESC LIMIT 1"
        ).fetchone()
        if ref is None:
            return {"reference": None, "stats": {}, "samples": [], "variants": []}

        ref_id = ref["id"]
        stats = {
            "valid": conn.execute(
                "SELECT COUNT(*) c FROM samples WHERE reference_id=? AND status='VALID'",
                (ref_id,),
            ).fetchone()["c"],
            "needs_reupload": conn.execute(
                "SELECT COUNT(*) c FROM samples WHERE reference_id=? "
                "AND status='NEEDS_REUPLOAD'",
                (ref_id,),
            ).fetchone()["c"],
            "stale": conn.execute(
                "SELECT COUNT(*) c FROM samples WHERE reference_id IS NOT NULL "
                "AND reference_id<>?",
                (ref_id,),
            ).fetchone()["c"],
        }
        samples = conn.execute(
            "SELECT s.*, "
            "  (SELECT COUNT(*) FROM sample_variants v WHERE v.sample_id=s.id) "
            "    AS variant_count, "
            "  (SELECT COUNT(*) FROM sample_variants v "
            "   WHERE v.sample_id=s.id AND v.gt <> '0/0') AS detected_count "
            "FROM samples s ORDER BY s.name"
        ).fetchall()
        # 首页汇总只统计检出记录（GT 0/0 表示已检测但未检出所列 ALT，
        # 不参与规范变异汇总与检出样本数）
        variants = conn.execute(
            "SELECT v.pos, v.ref, v.alt, v.kind, COUNT(DISTINCT v.sample_id) AS n_samples "
            "FROM sample_variants v "
            "JOIN samples s ON s.id = v.sample_id AND s.reference_id=? AND s.status='VALID' "
            "WHERE v.gt <> '0/0' "
            "GROUP BY v.pos, v.ref, v.alt "
            "ORDER BY v.pos, v.ref, v.alt",
            (ref_id,),
        ).fetchall()
        return {
            "reference": ref,
            "stats": stats,
            "samples": samples,
            "variants": variants,
        }


def get_variant_detail(pos: int, ref: str, alt: str) -> dict | None:
    with read_tx() as conn:
        current = conn.execute(
            "SELECT id FROM reference_sequences ORDER BY id DESC LIMIT 1"
        ).fetchone()
        if current is None:
            return None
        rows = conn.execute(
            "SELECT s.name AS sample, v.gt, v.orig_pos, v.orig_ref, v.orig_alt, v.kind "
            "FROM sample_variants v "
            "JOIN samples s ON s.id = v.sample_id "
            "WHERE v.pos=? AND v.ref=? AND v.alt=? "
            "AND v.gt <> '0/0' "
            "AND s.reference_id=? AND s.status='VALID' "
            "ORDER BY s.name",
            (pos, ref, alt, current["id"]),
        ).fetchall()
        if not rows:
            return None
        return {"pos": pos, "ref": ref, "alt": alt, "rows": rows}


# --------------------------------------------------------------------------- #
# 三联核查（双亲 + 子代，仅 SNP 位点）
# --------------------------------------------------------------------------- #
def run_trio_check(selection: dict[str, int]) -> dict:
    """在同一快照事务内完成一次三联核查。

    selection: {"father": sample_id, "mother": sample_id, "child": sample_id}
    三个样本必须互不相同，且都是当前最新参考下的 VALID 样本；
    所选样本已被重传（id 失效）或在参考重验后变为待重传时，抛出
    ValidationError 提示重新选择。

    整个核查在单个读事务（快照）内完成，保证「同一时刻、同一版参考下
    的有效数据」被一次性读取；上传/重验走写锁串行化，不会读到半更新状态。
    本函数只做只读核查，不持久化任何结果，因此样本重传或参考重验后再次
    提交即自然使用当前数据；失败的上传本身会整体回滚，原有核查依据保留。
    """
    ids = [selection.get(role) for role in ROLES]
    if any(sid is None for sid in ids):
        raise ValidationError("请为父亲、母亲、子代三个角色各选择一个样本。")
    if len(set(ids)) != 3:
        raise ValidationError("三联核查必须选择三个不同的样本，请重新选择。")

    with read_tx() as conn:
        current = conn.execute(
            "SELECT id, length, uploaded_at FROM reference_sequences "
            "ORDER BY id DESC LIMIT 1"
        ).fetchone()
        if current is None:
            raise ValidationError("尚未上传参考序列，无法核查。")
        ref_id = current["id"]

        chosen: dict[str, sqlite3.Row] = {}
        for role in ROLES:
            row = conn.execute(
                "SELECT id, name, status, reference_id FROM samples WHERE id=?",
                (selection[role],),
            ).fetchone()
            if row is None:
                raise ValidationError(
                    "所选样本已失效（可能已被删除或重传），请重新选择三个样本。"
                )
            if row["reference_id"] != ref_id or row["status"] != "VALID":
                raise ValidationError(
                    f"所选样本「{row['name']}」不是当前参考 #{ref_id} 下的有效样本"
                    "（参考重验后待重传或属于历史参考），请重新选择。"
                )
            chosen[role] = row

        # 仅取单碱基替换记录；0/0 行同样取出（它们是「已检测」的位点证据）
        records = conn.execute(
            "SELECT v.sample_id, s.name AS sample, v.pos, v.ref, v.alt, v.gt "
            "FROM sample_variants v "
            "JOIN samples s ON s.id = v.sample_id "
            "WHERE v.sample_id IN (?, ?, ?) "
            "AND length(v.ref)=1 AND length(v.alt)=1 AND v.ref <> v.alt "
            "ORDER BY v.pos",
            ids,
        ).fetchall()

    role_of_sample = {chosen[r]["id"]: r for r in ROLES}
    names = {r: chosen[r]["name"] for r in ROLES}
    flat = [
        {
            "sample": rec["sample"],
            "role": role_of_sample[rec["sample_id"]],
            "pos": rec["pos"],
            "ref": rec["ref"],
            "alt": rec["alt"],
            "gt": rec["gt"],
        }
        for rec in records
    ]
    result = summarize_trio_loci(flat, names)
    result.update(
        {
            "reference": dict(current),
            "names": names,
            "sample_ids": {r: chosen[r]["id"] for r in ROLES},
            "checked_at": datetime.now(timezone.utc)
            .strftime("%Y-%m-%d %H:%M:%S UTC"),
            "n_loci": len(result["consistent"])
            + len(result["violations"])
            + len(result["excluded"]),
        }
    )
    return result
