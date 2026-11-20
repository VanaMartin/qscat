"""Spin-separated CI overlaps with explicit fixed-core-hole sectors.

Canonical order is all alpha spatial orbitals then all beta spatial orbitals.
No native CSF/product conversion, Clebsch-Gordan phase inference or state solve.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import comb

import numpy as np


@dataclass(frozen=True)
class CoreSectorCI:
    """Root-by-alpha-string-by-beta-string coefficients in one core-hole sector."""

    ncore: int
    nactive: int
    core_alpha: tuple[int, ...]
    core_beta: tuple[int, ...]
    active_electrons: tuple[int, int]
    coefficients: np.ndarray


@lru_cache(maxsize=64)
def _active_strings(orbitals, electrons):
    strings = sorted(combinations(range(orbitals), electrons), key=lambda x: sum(1 << i for i in x))
    return tuple(strings)


def _validate(state):
    if not isinstance(state, CoreSectorCI):
        raise ValueError("Explicit CoreSectorCI required")
    if (
        type(state.ncore) is not int
        or not 0 <= state.ncore <= 4
        or type(state.nactive) is not int
        or not 0 <= state.nactive <= 11
        or len(state.active_electrons) != 2
    ):
        raise ValueError("Bounded core/active partition and two spin counts required")
    for core, active in zip(
        (state.core_alpha, state.core_beta), state.active_electrons, strict=True
    ):
        if (
            not isinstance(core, tuple)
            or any(type(i) is not int or not 0 <= i < state.ncore for i in core)
            or tuple(sorted(set(core))) != core
            or type(active) is not int
            or not 0 <= active <= state.nactive
            or len(core) + active > 8
        ):
            raise ValueError("Canonical core occupations and at most eight electrons/spin required")
    shape = tuple(comb(state.nactive, n) for n in state.active_electrons)
    c = np.asarray(state.coefficients)
    if (
        c.ndim != 3
        or c.shape[1:] != shape
        or not 1 <= len(c) <= 64
        or c.size > 20_000_000
        or not np.isfinite(c).all()
    ):
        raise ValueError("Finite bounded root-by-alpha-by-beta CI array required")
    return c


def _configuration(state):
    return state.core_alpha, state.core_beta, state.active_electrons


def _bra(coefficients):
    flat = coefficients.reshape(len(coefficients), -1)
    return flat.conj() if np.iscomplexobj(flat) else flat


class CoreSectorMetric:
    """Immutable spatial cross metric with bounded occupied-minor caches; no inverse."""

    def __init__(
        self, ncore, nactive, orbital_metric, *, max_pairs=1_000_000, cache_bytes=64 * 1024**2
    ):
        if (
            type(ncore) is not int
            or not 0 <= ncore <= 4
            or type(nactive) is not int
            or not 0 <= nactive <= 11
            or type(max_pairs) is not int
            or not 1 <= max_pairs <= 1_000_000
            or type(cache_bytes) is not int
            or not 1 <= cache_bytes <= 64 * 1024**2
        ):
            raise ValueError("Bounded partition, determinant-pair and cache caps required")
        matrix = np.asarray(orbital_metric)
        if matrix.shape != (ncore + nactive, ncore + nactive) or not np.isfinite(matrix).all():
            raise ValueError("Finite spatial core-plus-active cross metric required")
        self.ncore, self.nactive = ncore, nactive
        self.matrix = matrix.copy()
        self.matrix.setflags(write=False)
        self.max_pairs, self.cache_bytes = max_pairs, cache_bytes
        self._cache = OrderedDict()

    def _check(self, state):
        coefficients = _validate(state)
        if (state.ncore, state.nactive) != (self.ncore, self.nactive):
            raise ValueError("CI state and metric partitions differ")
        return coefficients

    def _minor(self, left_core, right_core, left_active, right_active):
        key = left_core, right_core, left_active, right_active
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        electrons = len(left_core) + left_active
        if electrons != len(right_core) + right_active:
            raise ValueError("Occupied minor spin-electron counts differ")
        left = [
            left_core + tuple(self.ncore + i for i in s)
            for s in _active_strings(self.nactive, left_active)
        ]
        right = [
            right_core + tuple(self.ncore + i for i in s)
            for s in _active_strings(self.nactive, right_active)
        ]
        pairs = len(left) * len(right)
        if pairs > self.max_pairs:
            raise ValueError("Occupied determinant pair cap exceeded")
        dtype = np.result_type(self.matrix.dtype, float)
        required = pairs * np.dtype(dtype).itemsize
        if required > self.cache_bytes:
            raise ValueError("One occupied-minor matrix exceeds cache cap")
        result = np.empty(pairs, dtype=dtype)
        if not electrons:
            result.fill(1)
        else:
            left, right = np.asarray(left, int), np.asarray(right, int)
            for start in range(0, pairs, 2048):
                indices = np.arange(start, min(pairs, start + 2048))
                a, b = left[indices // len(right)], right[indices % len(right)]
                result[start : start + len(indices)] = np.linalg.det(
                    self.matrix[a[:, :, None], b[:, None, :]]
                )
        result = result.reshape(len(left), len(right))
        if not np.isfinite(result).all():
            raise ValueError("Nonfinite occupied-spin minors")
        result.setflags(write=False)
        while (
            self._cache
            and sum(v.nbytes for v in self._cache.values()) + required > self.cache_bytes
        ):
            self._cache.popitem(last=False)
        self._cache[key] = result
        return result

    def _spin_metrics(self, left, right):
        totals = [
            tuple(
                len(c) + n
                for c, n in zip((s.core_alpha, s.core_beta), s.active_electrons, strict=True)
            )
            for s in (left, right)
        ]
        if totals[0] != totals[1]:
            return None
        return tuple(
            self._minor(a, b, na, nb)
            for a, b, na, nb in zip(
                (left.core_alpha, left.core_beta),
                (right.core_alpha, right.core_beta),
                left.active_electrons,
                right.active_electrons,
                strict=True,
            )
        )

    def overlaps(self, left: CoreSectorCI, right: CoreSectorCI) -> np.ndarray:
        """Root cross amplitudes, including all core-active occupied minors."""
        a, b = self._check(left), self._check(right)
        result = np.zeros((len(a), len(b)), dtype=np.result_type(a, b, self.matrix, float))
        kernels = self._spin_metrics(left, right)
        if kernels is None:
            return result
        alpha, beta = kernels
        for start in range(0, len(b), 8):
            transformed = np.stack([alpha @ v @ beta.T for v in b[start : start + 8]])
            result[:, start : start + len(transformed)] = (
                _bra(a) @ transformed.reshape(len(transformed), -1).T
            )
        return result


def remove_spin_orbital(state: CoreSectorCI, spin: int, orbital: int) -> CoreSectorCI | None:
    """Signed exterior-coordinate deletion; None is exact zero, not missing data."""
    coefficients = _validate(state)
    if (
        type(spin) is not int
        or spin not in (0, 1)
        or type(orbital) is not int
        or not 0 <= orbital < state.ncore + state.nactive
    ):
        raise ValueError("A valid spatial orbital and alpha0/beta1 spin required")
    cores = [state.core_alpha, state.core_beta]
    counts = list(state.active_electrons)
    alpha_total = len(state.core_alpha) + counts[0]
    preceding_spin = alpha_total if spin else 0
    if orbital < state.ncore:
        if orbital not in cores[spin]:
            return None
        position = cores[spin].index(orbital)
        cores[spin] = tuple(i for i in cores[spin] if i != orbital)
        result = (-1) ** (preceding_spin + position) * coefficients
    else:
        if counts[spin] == 0:
            return None
        local = orbital - state.ncore
        strings = _active_strings(state.nactive, counts[spin])
        counts[spin] -= 1
        destination = {s: i for i, s in enumerate(_active_strings(state.nactive, counts[spin]))}
        shape = (len(coefficients), *(comb(state.nactive, n) for n in counts))
        result = np.zeros(shape, dtype=coefficients.dtype)
        for source, occupied in enumerate(strings):
            if local not in occupied:
                continue
            position = occupied.index(local)
            target = destination[occupied[:position] + occupied[position + 1 :]]
            sign = (-1) ** (preceding_spin + len(cores[spin]) + position)
            if spin:
                result[:, :, target] = sign * coefficients[:, :, source]
            else:
                result[:, target, :] = sign * coefficients[:, source, :]
    return CoreSectorCI(state.ncore, state.nactive, cores[0], cores[1], tuple(counts), result)


def spin_lower(state: CoreSectorCI) -> CoreSectorCI | None:
    """Apply unnormalized S-minus to a paired frozen-core target, not native CG phases."""
    c = _validate(state)
    if state.core_alpha != tuple(range(state.ncore)) or state.core_beta != state.core_alpha:
        raise ValueError("Spin lowering currently requires a paired complete frozen core")
    na, nb = state.active_electrons
    if not na or nb == state.nactive:
        return None
    alpha, beta = _active_strings(state.nactive, na), _active_strings(state.nactive, nb)
    out_alpha = {s: i for i, s in enumerate(_active_strings(state.nactive, na - 1))}
    out_beta = {s: i for i, s in enumerate(_active_strings(state.nactive, nb + 1))}
    result = np.zeros((len(c), len(out_alpha), len(out_beta)), dtype=c.dtype)
    for i, a in enumerate(alpha):
        for position, orbital in enumerate(a):
            row = out_alpha[a[:position] + a[position + 1 :]]
            for j, b in enumerate(beta):
                if orbital in b:
                    continue
                insertion = sum(k < orbital for k in b)
                column = out_beta[tuple(sorted((*b, orbital)))]
                exponent = state.ncore + position + (state.ncore + na - 1) + state.ncore + insertion
                result[:, row, column] += (-1) ** exponent * c[:, i, j]
    return CoreSectorCI(
        state.ncore, state.nactive, state.core_alpha, state.core_beta, (na - 1, nb + 1), result
    )


def removal_overlaps(
    left: CoreSectorCI,
    right: CoreSectorCI,
    spin: int,
    metric: CoreSectorMetric,
    *,
    max_left_bytes=1024**3,
) -> np.ndarray:
    """Complete one-spin orbital/target removal Gram, grouped by core-hole sector.

    Left removal arrays are retained only within a checked byte budget; right
    arrays and transformed roots are processed one orbital/eight roots at a time.
    """
    a, b = metric._check(left), metric._check(right)
    if (
        type(spin) is not int
        or spin not in (0, 1)
        or type(max_left_bytes) is not int
        or not 1 <= max_left_bytes <= 1024**3
    ):
        raise ValueError("A spin and bounded left-removal memory budget required")
    counts = list(left.active_electrons)
    active_bytes = 0
    if counts[spin]:
        counts[spin] -= 1
        active_bytes = (
            left.nactive
            * len(a)
            * comb(left.nactive, counts[0])
            * comb(left.nactive, counts[1])
            * a.dtype.itemsize
        )
    core = left.core_alpha if spin == 0 else left.core_beta
    if active_bytes + len(core) * a.nbytes > max_left_bytes:
        raise ValueError("Left removal-sector memory budget exceeded")
    orbitals = left.ncore + left.nactive
    removed = [remove_spin_orbital(left, spin, i) for i in range(orbitals)]
    groups = {}
    for i, state in enumerate(removed):
        if state is not None:
            groups.setdefault(_configuration(state), []).append((i, state))
    result = np.zeros(
        (orbitals, len(a), orbitals, len(b)), dtype=np.result_type(a, b, metric.matrix, float)
    )
    for j in range(orbitals):
        right_state = remove_spin_orbital(right, spin, j)
        if right_state is None:
            continue
        for group in groups.values():
            kernels = metric._spin_metrics(group[0][1], right_state)
            if kernels is None:
                continue
            alpha, beta = kernels
            for start in range(0, len(b), 8):
                transformed = np.stack(
                    [alpha @ v @ beta.T for v in right_state.coefficients[start : start + 8]]
                )
                flat = transformed.reshape(len(transformed), -1).T
                for i, state in group:
                    result[i, :, j, start : start + len(transformed)] = (
                        _bra(state.coefficients) @ flat
                    )
    return result
