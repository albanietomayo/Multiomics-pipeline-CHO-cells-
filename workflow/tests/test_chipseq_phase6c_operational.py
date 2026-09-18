#!/usr/bin/env python3
"""Phase 6C operational safety tests; no scheduler or biological tools run."""

import importlib.util
import sys
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
SUBMIT_PATH = ROOT / "workflow/slurm/submit_chipseq_peak_calling.py"
RULES_PATH = ROOT / "workflow/rules/chipseq_peak_calling.smk"
SBATCH_PATH = ROOT / "workflow/slurm/chipseq_peak_calling.sbatch"
PPQT_ENV_PATH = ROOT / "workflow/envs/chipseq_phantompeakqualtools.yaml"

SPEC = importlib.util.spec_from_file_location("submit_peak_calling", SUBMIT_PATH)
SUBMIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUBMIT)


class Phase6COperationalTests(unittest.TestCase):
    def test_check_mode_cannot_reach_sbatch(self):
        with mock.patch.object(sys, "argv", [str(SUBMIT_PATH), "--check"]), \
                mock.patch.object(SUBMIT, "local_validation") as validation, \
                mock.patch.object(SUBMIT.subprocess, "check_output") as check_output, \
                mock.patch("builtins.print"):
            SUBMIT.main()
        validation.assert_called_once_with()
        check_output.assert_not_called()

    def test_failed_preflight_prevents_submit_side_effects(self):
        with mock.patch.object(sys, "argv", [str(SUBMIT_PATH), "--submit"]), \
                mock.patch.object(
                    SUBMIT, "local_validation", side_effect=ValueError("fail closed")
                ), \
                mock.patch.object(SUBMIT.tempfile, "mkdtemp") as mkdtemp, \
                mock.patch.object(SUBMIT.subprocess, "check_output") as check_output:
            with self.assertRaisesRegex(ValueError, "fail closed"):
                SUBMIT.main()
        mkdtemp.assert_not_called()
        check_output.assert_not_called()

    def test_submission_is_explicit_and_single_sbatch_site(self):
        source = SUBMIT_PATH.read_text(encoding="utf-8")
        self.assertIn('modes.add_argument("--check"', source)
        self.assertIn('modes.add_argument("--submit"', source)
        self.assertIn("mutually_exclusive_group(required=True)", source)
        self.assertEqual(source.count('"sbatch"'), 1)

    def test_ppqt_interface_and_atomic_rerun_contract(self):
        rules = RULES_PATH.read_text(encoding="utf-8")
        command = rules.split("rule chipseq_phantompeakqualtools:", 1)[1].split(
            "rule chipseq_parse_phantompeakqualtools:", 1
        )[0]
        for argument in ("-c=", "-p=", "-savp=", "-out="):
            self.assertIn(argument, command)
        self.assertNotIn("-clean", command)
        self.assertNotIn("-i=", command)
        self.assertIn("CHIP_PPQT_TABLE_TMP", command)
        self.assertIn("CHIP_PPQT_PLOT_TMP", command)
        self.assertIn('test -s "$CHIP_PPQT_TABLE_TMP"', command)
        self.assertIn('test -s "$CHIP_PPQT_PLOT_TMP"', command)

    def test_environment_pins_and_strict_channel_priority(self):
        environment = PPQT_ENV_PATH.read_text(encoding="utf-8")
        self.assertIn("phantompeakqualtools=1.2.2", environment)
        self.assertIn("samtools=1.24", environment)
        wrapper = SBATCH_PATH.read_text(encoding="utf-8")
        self.assertIn("export CONDA_CHANNEL_PRIORITY=strict", wrapper)


if __name__ == "__main__":
    unittest.main(verbosity=2)
