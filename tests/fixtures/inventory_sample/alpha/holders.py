"""Three holders that repeat the same save body."""

from dataclasses import dataclass


@dataclass(frozen=True)
class LevelsHolder:
    """Vibrational levels."""

    energies: list
    widths: list
    grid: str
    label: str

    def save(self, path):
        """Write the holder to path."""
        with open(path, "w") as fh:
            fh.write(repr(self))


@dataclass(frozen=True)
class StatesHolder:
    """Resonance states."""

    energies: list
    widths: list
    grid: str
    residuals: list

    def save(self, path):
        """Write the holder to path."""
        with open(path, "w") as fh:
            fh.write(repr(self))


@dataclass(frozen=True)
class CurveHolder:
    """A curve."""

    energies: list
    widths: list
    grid: str
    radii: list

    def save(self, path):
        """Write the holder to path."""
        with open(path, "w") as fh:
            fh.write(repr(self))
