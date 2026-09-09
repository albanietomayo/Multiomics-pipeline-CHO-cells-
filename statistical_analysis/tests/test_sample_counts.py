import sys
from pathlib import Path
import unittest
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from collect_sample_counts import collect, normalize
from plot_sample_counts import make_figure

class SampleCountsTests(unittest.TestCase):
    def fixture(self):
        return pd.DataFrame([("s1", " A ", "RNA-seq"), ("s1", "A", " RNA "),
                             ("s2", "A", "Multiome"), ("s2", "A", "ATAC"),
                             ("s3", "B", "scATAC"), ("s4", "B", "ChIP-seq")],
                            columns=["sample_id", "cell_line", "sequencing_type"])

    def test_dedup_multiome_grid_and_plot(self):
        inventory, counts, duplicates, unknown = collect(self.fixture())
        actual = {(r.cell_line, r.sequencing_type): r.count for r in counts.itertuples()}
        self.assertEqual(actual, {("A", "ATAC"): 1, ("A", "ChIP-seq"): 0, ("A", "RNA"): 2,
                                  ("B", "ATAC"): 1, ("B", "ChIP-seq"): 1, ("B", "RNA"): 0})
        self.assertEqual(duplicates, 2)
        self.assertEqual(unknown, ["ChIP-seq"])
        self.assertEqual(len(inventory), 5)
        self.assertEqual(counts["count"].sum(), len(inventory))
        fig = make_figure(counts)
        plotted = {(x, trace.name): int(y) for trace in fig.data for x, y in zip(trace.x, trace.y)}
        self.assertEqual(plotted, actual)

    def test_missing_rejected(self):
        for value in ("", " ", None, "NA"):
            frame = self.fixture()
            frame.loc[0, "cell_line"] = value
            with self.assertRaisesRegex(ValueError, "Missing required"):
                collect(frame)

    def test_conflict_rejected(self):
        frame = self.fixture()
        frame.loc[1, "cell_line"] = "B"
        with self.assertRaisesRegex(ValueError, "Conflicting.*s1"):
            collect(frame)

    def test_empty_rejected(self):
        with self.assertRaisesRegex(ValueError, "no sample"):
            collect(self.fixture().iloc[:0])

    def test_aliases(self):
        for label in ("RNA+ATAC", "atac / rna", "RNA-seq + ATAC-seq", "Multiome"):
            self.assertEqual(normalize(label), (["RNA", "ATAC"], False))
        self.assertEqual(normalize("scrna-seq"), (["RNA"], False))

if __name__ == "__main__":
    unittest.main()
