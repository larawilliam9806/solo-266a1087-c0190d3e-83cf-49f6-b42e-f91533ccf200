"""三联核查纯函数单元测试（无 Flask 依赖）。"""

import unittest

from app.trio import mendelian_compatible, summarize_trio_loci

NAMES = {"father": "爸爸", "mother": "妈妈", "child": "孩子"}


def rec(role, pos, ref, alt, gt):
    return {
        "sample": NAMES[role],
        "role": role,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "gt": gt,
    }


class MendelianCompatibleTest(unittest.TestCase):
    def test_basic_combinations(self):
        # 0/0 x 0/0 -> 只能 0/0
        self.assertTrue(mendelian_compatible("0/0", "0/0", "0/0"))
        self.assertFalse(mendelian_compatible("0/0", "0/0", "0/1"))
        self.assertFalse(mendelian_compatible("0/0", "0/0", "1/1"))
        # 0/0 x 1/1 -> 子代只能 0/1
        self.assertTrue(mendelian_compatible("0/0", "1/1", "0/1"))
        self.assertFalse(mendelian_compatible("0/0", "1/1", "0/0"))
        self.assertFalse(mendelian_compatible("0/0", "1/1", "1/1"))
        # 0/1 x 0/1 -> 0/0、0/1、1/1 均可
        for child in ("0/0", "0/1", "1/1"):
            self.assertTrue(mendelian_compatible("0/1", "0/1", child))
        # 1/1 x 1/1 -> 只能 1/1
        self.assertTrue(mendelian_compatible("1/1", "1/1", "1/1"))
        self.assertFalse(mendelian_compatible("1/1", "1/1", "0/1"))
        # 1/1 x 0/1 -> 0/1 或 1/1
        self.assertTrue(mendelian_compatible("1/1", "0/1", "0/1"))
        self.assertTrue(mendelian_compatible("1/1", "0/1", "1/1"))
        self.assertFalse(mendelian_compatible("1/1", "0/1", "0/0"))

    def test_zero_zero_treated_as_homozygous_ref(self):
        # 父 0/1 母 0/0：子代可为 0/0 或 0/1，不可能 1/1
        self.assertTrue(mendelian_compatible("0/1", "0/0", "0/0"))
        self.assertTrue(mendelian_compatible("0/1", "0/0", "0/1"))
        self.assertFalse(mendelian_compatible("0/1", "0/0", "1/1"))


class SummarizeTrioLociTest(unittest.TestCase):
    def test_violation_consistent_and_counts(self):
        records = [
            # pos1: 可解释（0/0 x 0/1 -> 0/0）
            rec("father", 1, "G", "C", "0/0"),
            rec("mother", 1, "G", "C", "0/1"),
            rec("child", 1, "G", "C", "0/0"),
            # pos2: 不可解释（0/0 x 0/0 生出 0/1）
            rec("father", 2, "A", "T", "0/0"),
            rec("mother", 2, "A", "T", "0/0"),
            rec("child", 2, "A", "T", "0/1"),
        ]
        out = summarize_trio_loci(records, NAMES)
        self.assertEqual([v["pos"] for v in out["consistent"]], [1])
        self.assertEqual([v["pos"] for v in out["violations"]], [2])
        self.assertEqual(out["excluded"], [])
        viol = out["violations"][0]
        self.assertEqual(viol["alts"], ["T"])
        self.assertEqual([r["role"] for r in viol["rows"]],
                         ["father", "mother", "child"])
        self.assertEqual(
            [r["gt"] for r in viol["rows"]], ["0/0", "0/0", "0/1"]
        )

    def test_missing_record_excluded_with_reason(self):
        records = [
            rec("father", 5, "C", "G", "0/1"),
            rec("child", 5, "C", "G", "1/1"),
            # 母亲缺该位点
        ]
        out = summarize_trio_loci(records, NAMES)
        self.assertEqual(out["consistent"], [])
        self.assertEqual(out["violations"], [])
        self.assertEqual(len(out["excluded"]), 1)
        item = out["excluded"][0]
        self.assertEqual(item["pos"], 5)
        self.assertIn("三方记录不齐全", item["reason"])
        self.assertIn(NAMES["mother"], item["reason"])
        mother_row = next(r for r in item["rows"] if r["role"] == "mother")
        self.assertIsNone(mother_row["gt"])

    def test_multiple_alts_excluded_even_if_mendel_otherwise_ok(self):
        records = [
            rec("father", 7, "T", "A", "0/1"),
            rec("mother", 7, "T", "C", "0/1"),       # 不同 ALT
            rec("child", 7, "T", "A", "0/0"),
        ]
        out = summarize_trio_loci(records, NAMES)
        self.assertEqual(out["violations"], [])
        self.assertEqual(len(out["excluded"]), 1)
        item = out["excluded"][0]
        self.assertEqual(item["alts"], ["A", "C"])
        self.assertIn("多等位基因", item["reason"])

    def test_zero_zero_alt_counts_toward_multi_alt(self):
        # 即使只是 0/0 行申报了不同 ALT，也算该位点出现多个 ALT
        records = [
            rec("father", 9, "G", "C", "1/1"),
            rec("mother", 9, "G", "A", "0/0"),
            rec("child", 9, "G", "C", "0/1"),
        ]
        out = summarize_trio_loci(records, NAMES)
        self.assertEqual(len(out["excluded"]), 1)
        self.assertEqual(out["excluded"][0]["alts"], ["A", "C"])

    def test_loci_sorted_and_empty_input(self):
        records = [
            rec("father", 12, "A", "C", "0/0"),
            rec("mother", 12, "A", "C", "0/0"),
            rec("child", 12, "A", "C", "0/0"),
            rec("father", 3, "A", "G", "0/0"),
            rec("mother", 3, "A", "G", "0/0"),
            rec("child", 3, "A", "G", "0/0"),
        ]
        out = summarize_trio_loci(records, NAMES)
        self.assertEqual([v["pos"] for v in out["consistent"]], [3, 12])
        empty = summarize_trio_loci([], NAMES)
        self.assertEqual(
            (empty["consistent"], empty["violations"], empty["excluded"]),
            ([], [], []),
        )


if __name__ == "__main__":
    unittest.main()
