"""Small offline tests for the notebooks' native ReGlyco adapter."""

import json
import re
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import reglyco_local as adapter

HELPER_DIR = Path(__file__).resolve().parent


def _notebook_code(name: str) -> str:
    """Return every code cell of a notebook as one Python source string."""

    document = json.loads((HELPER_DIR / name).read_text())
    return "\n".join(
        "".join(cell.get("source", []))
        for cell in document["cells"]
        if cell.get("cell_type") == "code"
    )


class LocalAdapterTests(unittest.TestCase):

    def test_release_manifest_pins_the_supported_binary(self) -> None:
        manifest = json.loads((HELPER_DIR / "reglyco_release.json").read_text())
        self.assertEqual(manifest["version"], "0.2.0")
        self.assertEqual(manifest["platform"], "linux-x86_64")
        self.assertEqual(manifest["asset"], "reglyco-0.2.0-linux-x86_64")
        self.assertIn("/v0.2.0/reglyco-0.2.0-linux-x86_64", manifest["url"])
        self.assertEqual(
            manifest["sha256"],
            "19b6d763deceae9f32b5dd1cc9c509d97075a250441122bff948e0af10df783b",
        )
        self.assertRegex(manifest["source_revision"], r"^[0-9a-f]{40}$")
        self.assertEqual(adapter.REGLYCO_VERSION, "0.2.0")

    def test_search_budget_args_default_to_auto(self) -> None:
        self.assertEqual(adapter.search_budget_args(), ["--search-budget", "auto"])
        self.assertEqual(adapter.search_budget_args("AUTO"), ["--search-budget", "auto"])
        self.assertEqual(adapter.search_budget_args("manual"), ["--search-budget", "manual"])
        with self.assertRaises(ValueError):
            adapter.search_budget_args("fast")

    def test_notebooks_use_the_auto_attachment_budget(self) -> None:
        for name in (
            "ReGlyco_Binder_Design_filter.ipynb",
            "Local_Conformational_Ensemble.ipynb",
        ):
            source = _notebook_code(name)
            self.assertIn("search_budget_args", source, name)
            self.assertTrue(re.search(r"search_budget_args\(", source), name)
            # The attachment search must never pin a manual budget again.
            self.assertNotIn("\"--population\"", source, name)
            self.assertNotIn("'--population'", source, name)
            self.assertNotIn("\"--generations\"", source, name)
            self.assertNotIn("'--generations'", source, name)
            self.assertNotIn("\"--search-budget\"", source, name)
            self.assertNotIn("'--search-budget'", source, name)

    def test_binder_notebook_keeps_the_glycoprotein_render_contract(self) -> None:
        source = _notebook_code("ReGlyco_Binder_Design_filter.ipynb")
        # --no-system suppresses the force-field bundle, not the PDB.
        self.assertIn("--no-system", source)
        self.assertIn("glycoprotein.pdb", source)
        self.assertIn("addGlycanHighlights", source)

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
