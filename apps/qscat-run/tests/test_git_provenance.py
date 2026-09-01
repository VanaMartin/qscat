"""The `git_sha` a manifest records is the binding between a published
artifact and the code that produced it -- the artifact store addresses
downloads by it, so a run that cannot determine it has produced an
unciteable result and must say so rather than writing `"unknown"` into a
plausible-looking manifest.

Three committed O2 sweeps (`validation/factory/results/o2-*-ve/`) shipped
`"git_sha": "unknown"` nine days after the baking mechanism landed, because
nothing here covered any of it.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import pytest
from qscat_run.artifacts import _git_sha

_HEX40 = re.compile(r"\A[0-9a-f]{40}\Z")

_DOCKERFILE = Path(__file__).resolve().parents[3] / "docker" / "Dockerfile"


def test_reads_the_repo_when_nothing_is_baked_in(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("QSCAT_GIT_SHA", raising=False)
    assert _HEX40.match(_git_sha())


def test_honours_a_baked_in_sha(monkeypatch: pytest.MonkeyPatch) -> None:
    """The containerised path: the build context excludes `.git`, so the host
    passes its SHA in as a build arg."""
    baked = "0" * 40
    monkeypatch.setenv("QSCAT_GIT_SHA", baked)
    assert _git_sha() == baked


@pytest.mark.parametrize("poison", ["unknown", "", "   ", "UNKNOWN", "not-a-sha", "HEAD"])
def test_a_non_sha_in_the_environment_does_not_defeat_the_fallback(
    monkeypatch: pytest.MonkeyPatch, poison: str
) -> None:
    """`ARG GIT_SHA=unknown` in the Dockerfile means an image built without
    `--build-arg GIT_SHA` carries the literal string `unknown` in the
    environment. Accepting any non-empty value made that sentinel outrank the
    working `git rev-parse` fallback, silently -- which is how the O2 sweeps
    recorded `"unknown"` while sitting in a readable repo."""
    monkeypatch.setenv("QSCAT_GIT_SHA", poison)
    assert _HEX40.match(_git_sha())


def test_warns_but_proceeds_when_the_sha_is_genuinely_undeterminable(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """No baked SHA and no repo to read.

    This used to raise, because published artifacts were addressed by commit
    and a run without one could not be published at all. They are addressed by
    content digest now, so a missing SHA costs traceability and breaks
    nothing -- a local figure should not fail over forty bytes of metadata.
    The caller still has to hear about it, or this regresses to the silence
    that let three sweeps ship `"unknown"` unnoticed.
    """
    monkeypatch.setenv("QSCAT_GIT_SHA", "unknown")
    monkeypatch.delenv("QSCAT_ALLOW_UNKNOWN_SHA", raising=False)
    monkeypatch.setattr("qscat_run.artifacts._REPO_PROBE_DIR", tmp_path)
    with pytest.warns(RuntimeWarning, match="cannot determine the commit SHA"):
        assert _git_sha() == "unknown"


def test_the_warning_can_be_silenced_deliberately(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Running from an unpacked tarball is legitimate; the opt-out says so
    once rather than warning on every run."""
    monkeypatch.setenv("QSCAT_GIT_SHA", "unknown")
    monkeypatch.setenv("QSCAT_ALLOW_UNKNOWN_SHA", "1")
    monkeypatch.setattr("qscat_run.artifacts._REPO_PROBE_DIR", tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _git_sha() == "unknown"


def test_every_dockerfile_stage_that_can_run_qscat_run_stamps_the_sha() -> None:
    """`ENV` does not cross a `FROM` boundary and `COPY --from=` copies files,
    not environment. `runtime` starts from a fresh upstream base, so it needs
    its own `ARG`/`ENV` pair -- measured absent in a built `qmodeling:runtime`
    while `qmodeling:test-deps` had it."""
    text = _DOCKERFILE.read_text()
    stages = {
        name: body
        for name, body in zip(
            re.findall(r"^FROM .* AS (\S+)", text, re.M),
            re.split(r"^FROM .* AS \S+", text, flags=re.M)[1:],
            strict=True,
        )
    }
    fresh_base = [
        name
        for name in stages
        if re.search(rf"^FROM (?!build|test-deps|test\b).*AS {name}$", text, re.M)
    ]
    assert "runtime" in fresh_base, "guard assumes runtime starts from a fresh base"
    for name in fresh_base:
        assert "QSCAT_GIT_SHA" in stages[name], (
            f"stage {name!r} starts from a fresh base, so it inherits no ENV; "
            "it needs its own ARG GIT_SHA / ENV QSCAT_GIT_SHA"
        )
