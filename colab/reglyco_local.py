"""Small local adapter for the released ReGlyco executable.

The notebooks deliberately keep provider level fixed internally.  Users choose
structures, sites, seeds, and workflow settings; they do not select an asset
level.  The helper is also usable with a locally built binary while a release
is being staged.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


_RELEASE_METADATA_PATH = Path(__file__).with_name("reglyco_release.json")


def _release_metadata() -> dict[str, Any]:
    try:
        value = json.loads(_RELEASE_METADATA_PATH.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


_RELEASE = _release_metadata()
REGLYCO_VERSION = os.environ.get("REGLYCO_VERSION", str(_RELEASE.get("version", "0.2.0")))
REGLYCO_API_BASE = os.environ.get("REGLYCO_API_BASE", "https://glycoshape.io")
_FIXED_PROVIDER_LEVEL = "1"


class ReGlycoCommandError(RuntimeError):
    """A ReGlyco command failed before producing a scientific result."""

    def __init__(self, message: str, result: "CommandResult") -> None:
        super().__init__(message)
        self.result = result


@dataclass(frozen=True)
class CommandResult:
    command: tuple[str, ...]
    output_dir: Path
    returncode: int
    stdout: str
    stderr: str

    @property
    def ok(self) -> bool:
        return self.returncode == 0

    @property
    def command_text(self) -> str:
        return " ".join(_quote(arg) for arg in self.command)


def _quote(value: str) -> str:
    if value and all(ch.isalnum() or ch in ".-_/:=@" for ch in value):
        return value
    return repr(value)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=120) as response, temporary.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def _build_from_checkout(root: Path, destination: Path) -> Path:
    manifest = root / "Cargo.toml"
    if not manifest.is_file():
        raise FileNotFoundError(f"ReGlyco checkout has no Cargo.toml: {root}")
    subprocess.run(
        ["cargo", "build", "--release", "--locked", "--manifest-path", str(manifest)],
        check=True,
        cwd=root,
    )
    built = root / "target" / "release" / "reglyco"
    if not built.is_file():
        raise FileNotFoundError(f"Cargo build did not produce {built}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built, destination)
    destination.chmod(0o755)
    return destination


def ensure_reglyco_binary(cache_dir: str | Path = ".reglyco-colab") -> Path:
    """Resolve a pinned release binary, local override, or source fallback."""

    override = os.environ.get("REGLYCO_BIN")
    if override:
        path = Path(override).expanduser().resolve()
        if not path.is_file() or not os.access(path, os.X_OK):
            raise FileNotFoundError(f"REGLYCO_BIN is not executable: {path}")
        return path

    cache = Path(cache_dir).expanduser().resolve()
    binary = cache / f"reglyco-{REGLYCO_VERSION}-linux-x86_64"
    expected = os.environ.get("REGLYCO_SHA256", str(_RELEASE.get("sha256", ""))).strip().lower()
    if binary.is_file() and (not expected or _sha256(binary) == expected):
        binary.chmod(0o755)
        return binary

    explicit_release_url = os.environ.get("REGLYCO_BINARY_URL", "").strip()
    release_url = explicit_release_url or str(_RELEASE.get("url", "")).strip()
    if release_url:
        try:
            _download(release_url, binary)
            binary.chmod(0o755)
            actual = _sha256(binary)
            if expected and actual != expected:
                binary.unlink(missing_ok=True)
                raise RuntimeError(f"ReGlyco checksum mismatch: expected {expected}, got {actual}")
            return binary
        except Exception:
            # A checked-in release manifest is intentionally usable before the
            # GitHub asset is uploaded.  Explicit URLs remain strict; the
            # automatic default may fall back to a locked Cargo install.
            if explicit_release_url:
                raise

    checkout = os.environ.get("REGLYCO_SOURCE_ROOT", "").strip()
    if checkout:
        return _build_from_checkout(Path(checkout).expanduser().resolve(), binary)

    cargo = shutil.which("cargo")
    if cargo is None:
        raise RuntimeError("ReGlyco is unavailable: set REGLYCO_BIN or install Rust/Cargo.")
    install_root = cache / "cargo-install"
    install_root.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            cargo,
            "install",
            "reglyco",
            "--locked",
            "--version",
            REGLYCO_VERSION,
            "--root",
            str(install_root),
        ],
        check=True,
    )
    binary = install_root / "bin" / "reglyco"
    if not binary.is_file():
        raise FileNotFoundError(f"cargo install did not produce {binary}")
    binary.chmod(0o755)
    return binary


def _common_provider_args(cache_dir: Path) -> list[str]:
    return [
        "--api-base",
        REGLYCO_API_BASE,
        "--cache",
        str(cache_dir / "assets"),
        "--anomer",
        "beta",
        "--level",
        _FIXED_PROVIDER_LEVEL,
    ]


def search_budget_args(mode: str = "auto") -> list[str]:
    """Return the flags selecting a ReGlyco attachment-search budget.

    `build` and `ensemble` accept `--search-budget auto|manual`.  Auto derives
    the population and generation count from the loaded site/conformer conflict
    graph, so the notebooks must not hard-code `--population`/`--generations`
    for them.  (`scan` keeps its own fast 32x25 Cookbook defaults.)
    """

    normalized = str(mode).strip().lower()
    if normalized not in {"auto", "manual"}:
        raise ValueError(f"Unsupported ReGlyco search budget mode: {mode!r}")
    return ["--search-budget", normalized]


def run_reglyco(
    arguments: Sequence[str],
    output_dir: str | Path,
    *,
    cache_dir: str | Path = ".reglyco-colab",
    seed: int | None = None,
    threads: int | None = 1,
    check: bool = False,
) -> CommandResult:
    """Run one local command with fixed provider settings and provenance."""

    output = Path(output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    cache = Path(cache_dir).expanduser().resolve()
    binary = ensure_reglyco_binary(cache)
    command = [str(binary), *arguments]
    if seed is not None and "--seed" not in command:
        command.extend(["--seed", str(int(seed))])
    # `scan` deliberately has no worker-thread option in the native CLI.  It
    # is a small deterministic topology/accessibility pass, while `build` and
    # `ensemble` expose the threaded prepared-scoring paths.  Do not append a
    # shared adapter option to a subcommand that does not accept it: clap
    # rejects the whole command with exit code 2 before producing scan.json.
    subcommand = next((arg for arg in command[1:] if not arg.startswith("-")), "")
    if threads is not None and subcommand in {"build", "ensemble"} and "--threads" not in command:
        command.extend(["--threads", str(max(1, int(threads)))])
    command.extend(_common_provider_args(cache))
    completed = subprocess.run(
        command,
        cwd=output,
        text=True,
        capture_output=True,
        check=False,
    )
    result = CommandResult(tuple(command), output, completed.returncode, completed.stdout, completed.stderr)
    (output / "reglyco.stdout.log").write_text(result.stdout)
    (output / "reglyco.stderr.log").write_text(result.stderr)
    (output / "provenance.json").write_text(
        json.dumps(provenance(result, seed=seed), indent=2) + "\n"
    )
    if check and not result.ok:
        raise ReGlycoCommandError(
            f"ReGlyco command failed ({result.returncode}): {result.command_text}", result
        )
    return result


def read_json(output_dir: str | Path, filename: str) -> dict[str, Any] | None:
    path = Path(output_dir) / filename
    if not path.is_file():
        return None
    try:
        value = json.loads(path.read_text())
    except json.JSONDecodeError:
        return None
    return value if isinstance(value, dict) else None


def read_report(output_dir: str | Path) -> dict[str, Any] | None:
    return read_json(output_dir, "report.json") or read_json(output_dir, "search.json")


def clash_status(output_dir: str | Path) -> str | None:
    report = read_report(output_dir) or {}
    status = report.get("clash_status") or report.get("clashStatus")
    if isinstance(status, dict):
        status = status.get("status")
    return str(status) if status is not None else None


def is_clash_free(output_dir: str | Path) -> bool:
    return clash_status(output_dir) in {"clash_free", "ClashFree"}


def output_structure(output_dir: str | Path) -> Path | None:
    for name in ("glycoprotein.pdb", "structure.pdb", "output.pdb", "ensemble.pdb"):
        path = Path(output_dir) / name
        if path.is_file() and path.stat().st_size > 0:
            return path
    return None


def provenance(result: CommandResult, *, seed: int | None = None) -> dict[str, Any]:
    return {
        "executable": str(result.command[0]),
        "version": REGLYCO_VERSION,
        "sha256": _sha256(Path(result.command[0])) if Path(result.command[0]).is_file() else None,
        "command": list(result.command),
        "command_text": result.command_text,
        "returncode": result.returncode,
        "seed": seed,
        "provider": {"api_base": REGLYCO_API_BASE, "level": _FIXED_PROVIDER_LEVEL, "anomer": "beta"},
    }


def assert_linux_x86_64() -> None:
    if platform.system() != "Linux" or platform.machine().lower() not in {"x86_64", "amd64"}:
        raise RuntimeError("The pinned Colab executable is built for Linux x86-64.")


__all__ = [
    "CommandResult",
    "ReGlycoCommandError",
    "assert_linux_x86_64",
    "clash_status",
    "ensure_reglyco_binary",
    "is_clash_free",
    "output_structure",
    "provenance",
    "read_json",
    "read_report",
    "run_reglyco",
    "search_budget_args",
]
