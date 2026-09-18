"""Small offline tests for the notebooks' native ReGlyco adapter."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import reglyco_local as adapter


class LocalAdapterTests(unittest.TestCase):
    def test_provider_level_and_seed_are_added_once(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary = root / "reglyco"
            binary.write_bytes(b"test executable")
            output = root / "output"
            captured = {}

            def fake_run(command, **kwargs):
                captured["command"] = list(command)
                return adapter.subprocess.CompletedProcess(command, 0, "ok\n", "")

            with patch.object(adapter, "ensure_reglyco_binary", return_value=binary), patch.object(
                adapter.subprocess, "run", side_effect=fake_run
            ):
                result = adapter.run_reglyco(
                    ["build", "--output", str(output)],
                    output,
                    cache_dir=root / "cache",
                    seed=0,
                    threads=3,
                )

            command = captured["command"]
            self.assertTrue(result.ok)
            self.assertEqual(command.count("--level"), 1)
            self.assertEqual(command[command.index("--level") + 1], "1")
            self.assertEqual(command[command.index("--seed") + 1], "0")
            self.assertEqual(command[command.index("--threads") + 1], "3")
            provenance = json.loads((output / "provenance.json").read_text())
            self.assertEqual(provenance["seed"], 0)
            self.assertEqual(provenance["provider"]["level"], "1")

    def test_existing_seed_and_threads_are_not_duplicated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary = root / "reglyco"
            binary.write_bytes(b"test executable")
            output = root / "output"
            captured = {}

            def fake_run(command, **kwargs):
                captured["command"] = list(command)
                return adapter.subprocess.CompletedProcess(command, 0, "", "")

            with patch.object(adapter, "ensure_reglyco_binary", return_value=binary), patch.object(
                adapter.subprocess, "run", side_effect=fake_run
            ):
                adapter.run_reglyco(
                    ["ensemble", "--seed", "12", "--threads", "5", "--output", str(output)],
                    output,
                    cache_dir=root / "cache",
                    seed=99,
                    threads=7,
                )

            command = captured["command"]
            self.assertEqual(command.count("--seed"), 1)
            self.assertEqual(command[command.index("--seed") + 1], "12")
            self.assertEqual(command.count("--threads"), 1)
            self.assertEqual(command[command.index("--threads") + 1], "5")

    def test_scan_does_not_receive_unsupported_threads_option(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary = root / "reglyco"
            binary.write_bytes(b"test executable")
            output = root / "output"
            captured = {}

            def fake_run(command, **kwargs):
                captured["command"] = list(command)
                return adapter.subprocess.CompletedProcess(command, 0, "", "")

            with patch.object(adapter, "ensure_reglyco_binary", return_value=binary), patch.object(
                adapter.subprocess, "run", side_effect=fake_run
            ):
                adapter.run_reglyco(
                    ["scan", "--protein", "input.pdb", "--output", str(output)],
                    output,
                    cache_dir=root / "cache",
                    seed=0,
                    threads=8,
                )

            command = captured["command"]
            self.assertEqual(command[1], "scan")
            self.assertNotIn("--threads", command)
            self.assertEqual(command[command.index("--level") + 1], "1")

    def test_report_and_output_classification(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "search.json").write_text(json.dumps({"clash_status": "clash_free"}))
            (root / "glycoprotein.pdb").write_text("ATOM\n")
            self.assertTrue(adapter.is_clash_free(root))
            self.assertEqual(adapter.clash_status(root), "clash_free")
            self.assertEqual(adapter.output_structure(root), root / "glycoprotein.pdb")


if __name__ == "__main__":
    unittest.main()
