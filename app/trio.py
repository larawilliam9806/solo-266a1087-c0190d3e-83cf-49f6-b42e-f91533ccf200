"""三联（双亲 + 子代）单碱基位点核查的纯函数，不依赖 Flask / SQLite。

核查范围仅为单碱基替换（SNP）位点。每份样本在同一规范 SNP 位点最多一条
记录（重复记录在上传阶段已被整份拒绝），记录的 GT 含义：

  * 0/0：该位点已检测、未检出所申报 ALT（纯合参考，两个等位基因均为 0）；
  * 0/1：杂合，一个参考等位基因 0 + 一个 ALT 等位基因 1；
  * 1/1：纯合 ALT，两个等位基因均为 1。

对「三方记录齐全、且三方在该位点只出现一种 ALT」的位点，依据孟德尔遗传
检验：子代的两个等位基因必须能分别从父亲、母亲各取得一个。其余位点
（任何一方缺记录，或三方记录中出现多个不同 ALT）单独列出并说明原因，
不做遗传判断。
"""

from __future__ import annotations

from typing import Iterable

# 角色键 -> 中文展示名
ROLE_LABELS = {"father": "父亲", "mother": "母亲", "child": "子代"}
ROLES = ("father", "mother", "child")

# 每条记录预期形如：
# {sample: 样本名, role: father/mother/child, pos, ref, alt, gt}


def gt_alleles(gt: str) -> tuple[str, str]:
    """把 GT 文本拆成两个等位基因（'0' 参考 / '1' ALT）。"""
    a, b = gt.split("/", 1)
    return a, b


def mendelian_compatible(
    father_gt: str, mother_gt: str, child_gt: str
) -> bool:
    """判断子代 GT 能否由双亲各传一个等位基因组合得到。

    双亲的四个单倍型（父亲两等位基因 × 母亲两等位基因）构成所有可能的
    子代无序基因型，只要子代 GT 在其中即视为可解释。
    """
    fa, fb = gt_alleles(father_gt)
    ma, mb = gt_alleles(mother_gt)
    ca, cb = gt_alleles(child_gt)
    child = frozenset((ca, cb))
    for f in (fa, fb):
        for m in (ma, mb):
            if frozenset((f, m)) == child:
                return True
    return False


def summarize_trio_loci(
    records: Iterable[dict], names: dict[str, str]
) -> dict:
    """把三方 SNP 记录按规范位点汇总并分类。

    参数
      records: 三方在当前参考下的全部 SNP 记录，每条含
               sample/role/pos/ref/alt/gt；
      names: {role: 样本名}，用于输出与缺记录提示。

    返回
      {
        "consistent": [位点...],   # 三方齐全、单一 ALT、遗传可解释
        "violations": [位点...],   # 三方齐全、单一 ALT、遗传不可解释
        "excluded":   [位点...],   # 缺记录或多个 ALT，附 reason
      }
    每个位点条目含 pos/ref/alts（有序去重 ALT 列表）、
    gts/samples（按 father/mother/child 排序）。
    """
    loci: dict[int, dict] = {}
    for rec in records:
        locus = loci.setdefault(
            rec["pos"],
            {"pos": rec["pos"], "ref": rec["ref"], "by_role": {}},
        )
        # 同一样本同一位点不可能有多条（上传时重复校验已拦截）
        locus["by_role"][rec["role"]] = rec

    consistent: list[dict] = []
    violations: list[dict] = []
    excluded: list[dict] = []

    for pos in sorted(loci):
        locus = loci[pos]
        by_role = locus["by_role"]
        present_roles = [r for r in ROLES if r in by_role]

        if len(present_roles) < 3:
            missing_roles = [r for r in ROLES if r not in by_role]
            missing_desc = "、".join(
                f"{ROLE_LABELS[r]}「{names[r]}」" for r in missing_roles
            )
            excluded.append(
                {
                    **_entry(locus, by_role),
                    "reason": f"三方记录不齐全：{missing_desc} 缺少该位点记录。",
                }
            )
            continue

        alts = sorted({by_role[r]["alt"] for r in ROLES})
        entry = _entry(locus, by_role, alts=alts)
        if len(alts) > 1:
            entry["reason"] = (
                f"三方在该位点出现 {len(alts)} 种不同 ALT"
                f"（{ '、'.join(alts) }），存在多等位基因情况，"
                "无法按单一 ALT 做亲子传递判断。"
            )
            excluded.append(entry)
            continue

        father_gt = by_role["father"]["gt"]
        mother_gt = by_role["mother"]["gt"]
        child_gt = by_role["child"]["gt"]
        if mendelian_compatible(father_gt, mother_gt, child_gt):
            consistent.append(entry)
        else:
            violations.append(entry)

    return {
        "consistent": consistent,
        "violations": violations,
        "excluded": excluded,
    }


def _entry(locus: dict, by_role: dict, alts: list[str] | None = None) -> dict:
    """构造模板可直接渲染的位点条目。"""
    if alts is None:
        alts = sorted({r["alt"] for r in by_role.values()})
    rows = []
    for role in ROLES:
        rec = by_role.get(role)
        rows.append(
            {
                "role": role,
                "role_label": ROLE_LABELS[role],
                "sample": rec["sample"] if rec else None,
                "gt": rec["gt"] if rec else None,
                "alt": rec["alt"] if rec else None,
            }
        )
    return {
        "pos": locus["pos"],
        "ref": locus["ref"],
        "alts": alts,
        "rows": rows,
    }
