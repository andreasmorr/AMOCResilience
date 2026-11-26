"""Check CLIMBER-X synthesis provenance, censoring, and displayed units."""
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

    def test_all_levels_bounds_and_independent_local_fit(self):
        raw = pd.read_csv(self.csv)
        frame = synthesis.load_climberx()
        self.assertEqual(len(frame), 135)
        self.assertEqual(set(frame.co2_ppm), set(synthesis.CLIMBERX_LEVELS))
        mask = frame.measure.eq("basin_volume")
        np.testing.assert_allclose(frame.loc[mask, ["value", "lower", "upper"]],
                                   raw.loc[mask, ["value", "lower", "upper"]] * 16)
        self.assertEqual(frame.provisional.sum(), 16)
        self.assertFalse(frame.loc[frame.measure.eq("local_resilience"), "provisional"].any())
        self.assertTrue(frame.loc[(frame.co2_ppm == 360) & frame.measure.eq("min_critical_shock"), "value"].isna().all())
        self.assertEqual(frame.attrs["analysis"]["n_short_grid"], 4)

    def test_modified_csv_rejected_by_provenance_hash(self):
        self.csv.write_text(self.csv.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            synthesis.load_climberx()

    def test_false_certainty_rejected_even_with_matching_hash(self):
        frame = pd.read_csv(self.csv)
        frame.loc[frame.co2_ppm.eq(400) & frame.measure.eq("basin_volume"), "provisional"] = False
        self.rewrite(frame)
        with self.assertRaisesRegex(ValueError, "provisional flags"):
            synthesis.load_climberx()

    def test_wrong_bounds_rejected_even_with_matching_hash(self):
        frame = pd.read_csv(self.csv)
        mask = frame.co2_ppm.eq(570) & frame.measure.eq("basin_volume")
        frame.loc[mask, "upper"] = frame.loc[mask, "lower"]
        self.rewrite(frame)
        with self.assertRaisesRegex(ValueError, "bounds disagree"):
            synthesis.load_climberx()

    def test_provisional_hollow_markers_and_area_errorbars(self):
        frame = synthesis.load_climberx()
        fig, axis = plt.subplots()
        self.addCleanup(plt.close, fig)
        synthesis._plot_panel(axis, None, None, "basin_volume", None, "area", "area",
                              None, None, None, frame)
        hollow = [line for line in axis.lines
                  if line.get_markerfacecolor() == "white" and len(line.get_xdata())]
        self.assertEqual(set(hollow[0].get_xdata()), {360., 400., 570., 600.})
        bars = [segment for collection in axis.collections for segment in collection.get_segments()]
        self.assertEqual(len(bars), 4)
        np.testing.assert_allclose([segment[1, 1] - segment[0, 1] for segment in bars], 16 / 81)


if __name__ == "__main__":
    unittest.main()
