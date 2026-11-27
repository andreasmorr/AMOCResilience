"""Check CLIMBER-X synthesis provenance, branch labels, and displayed units."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
UMBRELLA = Path(os.environ.get("AMOC_UMBRELLA", ROOT))
with patch.object(sys, "argv", ["synthesis_figure.py", "--umbrella", str(UMBRELLA)]):
    spec = importlib.util.spec_from_file_location("synthesis_under_test", ROOT / "synthesis_figure.py")
    synthesis = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(synthesis)


class ClimberSynthesisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.bundle = Path(self.temp.name)
        original = synthesis.CLIMBERX_CSV
        for name in (original.name, "analysis_metadata.json", "run_returns.csv"):
            shutil.copy2(original.parent / name, self.bundle / name)
        self.csv = self.bundle / original.name
        self.patch = patch.object(synthesis, "CLIMBERX_CSV", self.csv)
        self.patch.start()
        self.addCleanup(self.patch.stop)

    def rewrite(self, frame):
        frame.to_csv(self.csv, index=False)
        path = self.bundle / "analysis_metadata.json"
        metadata = json.loads(path.read_text())
        metadata["output_sha256"][self.csv.name] = hashlib.sha256(self.csv.read_bytes()).hexdigest()
        path.write_text(json.dumps(metadata))

    def test_all_levels_and_area_units(self):
        raw = pd.read_csv(self.csv)
        frame = synthesis.load_climberx()
        self.assertEqual(len(frame), 135)
        self.assertEqual(set(frame.co2_ppm), set(synthesis.CLIMBERX_LEVELS))
        mask = frame.measure.eq("basin_volume")
        np.testing.assert_allclose(frame.loc[mask, "value"], raw.loc[mask, "value"] * 16)
        local = frame.measure.eq("local_resilience")
        np.testing.assert_allclose(frame.loc[local, "value"], raw.loc[local, "value"])
        self.assertTrue(frame.loc[(frame.co2_ppm == 360) & frame.measure.eq("min_critical_shock"), "value"].isna().all())
        self.assertEqual(frame.attrs["analysis"]["n_grid"], 1215)

    def test_modified_csv_rejected_by_provenance_hash(self):
        self.csv.write_text(self.csv.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            synthesis.load_climberx()

    def test_wrong_branch_label_rejected_even_with_matching_hash(self):
        frame = pd.read_csv(self.csv)
        frame.loc[frame.co2_ppm.eq(330) & frame.measure.eq("basin_volume"), "attractor"] = "strong"
        self.rewrite(frame)
        with self.assertRaisesRegex(ValueError, "branch labels"):
            synthesis.load_climberx()

    def test_filled_markers_without_errorbars(self):
        frame = synthesis.load_climberx()
        fig, axis = plt.subplots()
        self.addCleanup(plt.close, fig)
        synthesis._plot_panel(axis, None, None, "basin_volume", None, "area", "area",
                              None, None, None, frame)
        self.assertFalse([line for line in axis.lines if line.get_markerfacecolor() == "white"])
        self.assertFalse([segment for collection in axis.collections
                          for segment in collection.get_segments()])
        plotted = sorted(x for line in axis.lines for x in line.get_xdata())
        self.assertEqual(plotted, sorted(float(co2) for co2 in synthesis.CLIMBERX_LEVELS))


if __name__ == "__main__":
    unittest.main()
