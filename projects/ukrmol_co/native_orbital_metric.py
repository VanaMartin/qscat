"""Loader-aware one-body Gaussian self metrics, not a native state comparator.

See docs/physics/co-native-orbital-conventions.md. All amplitudes are explicit;
the original export and MO coefficients are never modified or whitened.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import islice
from math import comb, factorial, gamma, pi

import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from scipy.special import gammainc, gammaincc, sph_harm_y


@dataclass(frozen=True)
class SphericalGaussian:
    """Native real-solid-harmonic shell; amplitudes already contain every norm."""

    center: tuple[float, float, float]
    angular: int
    exponents: tuple[float, ...]
    amplitudes: tuple[complex, ...]
    continuum: bool = False


def _shell(shell):
    if not isinstance(shell, SphericalGaussian):
        raise ValueError("Explicit SphericalGaussian shell required")
    center, exponents, amplitudes = map(
        np.asarray, [shell.center, shell.exponents, shell.amplitudes]
    )
    if (
        center.shape != (3,)
        or not np.isrealobj(center)
        or not np.isfinite(center).all()
        or type(shell.angular) is not int
        or not 0 <= shell.angular <= 10
        or exponents.ndim != 1
        or not 1 <= len(exponents) <= 128
        or not np.isrealobj(exponents)
        or not np.isfinite(exponents).all()
        or np.any(exponents <= 0)
        or amplitudes.shape != exponents.shape
        or not np.isfinite(amplitudes).all()
        or type(shell.continuum) is not bool
    ):
        raise ValueError("Finite bounded spherical Gaussian shell required")
    if shell.continuum and np.any(center != 0):
        raise ValueError("Native continuum shell must be centered at the origin")


def _basis(shells):
    shells = tuple(islice(iter(shells), 129))
    if not 1 <= len(shells) <= 128:
        raise ValueError("One to128 shells required")
    for shell in shells:
        _shell(shell)
    sizes = [2 * shell.angular + 1 for shell in shells]
    if sum(sizes) > 512 or sum(len(s.exponents) for s in shells) > 256:
        raise ValueError("AO or primitive basis cap exceeded")
    return shells, np.r_[0, np.cumsum(sizes)]


def solid_harmonics(angular: int, points: np.ndarray) -> np.ndarray:
    """R_lm in m=-l,...,+l order, via normalized complex harmonics (with CS phase)."""
    points = np.asarray(points)
    if type(angular) is not int or not 0 <= angular <= 10:
        raise ValueError("Angular degree must be an integer from0 to10")
    if points.ndim != 2 or points.shape[1] != 3 or len(points) > 100_000:
        raise ValueError("At most100000 three-dimensional points required")
    if not np.isrealobj(points) or not np.isfinite(points).all():
        raise ValueError("Finite real points required")
    radius = np.linalg.norm(points, axis=1)
    cosine = np.divide(points[:, 2], radius, out=np.zeros_like(radius), where=radius != 0)
    theta, phi = np.arccos(np.clip(cosine, -1, 1)), np.arctan2(points[:, 1], points[:, 0])
    values = np.empty((len(points), 2 * angular + 1))
    scale = np.sqrt(4 * pi / (2 * angular + 1)) * radius**angular
    values[:, angular] = scale * sph_harm_y(angular, 0, theta, phi).real
    for m in range(1, angular + 1):
        harmonic = np.sqrt(2) * (-1) ** m * scale * sph_harm_y(angular, m, theta, phi)
        values[:, angular + m], values[:, angular - m] = harmonic.real, harmonic.imag
    if not np.isfinite(values).all():
        raise ValueError("Unsupported solid-harmonic coordinate scale")
    return values


@lru_cache(maxsize=11)
def _polynomials(angular):
    """Independent Cartesian coefficients from the explicit Legendre polynomial."""
    powers = np.array(
        [(x, y, angular - x - y) for x in range(angular + 1) for y in range(angular - x + 1)]
    )
    lookup = {tuple(power): i for i, power in enumerate(powers)}
    matrix = np.zeros((2 * angular + 1, len(powers)))
    for m in range(angular + 1):
        polynomial = np.zeros(len(powers), complex)
        scale = 1 if m == 0 else np.sqrt(2 * factorial(angular - m) / factorial(angular + m))
        for k in range((angular - m) // 2 + 1):
            legendre = (
                (-1) ** k
                * factorial(2 * angular - 2 * k)
                / (
                    2**angular
                    * factorial(k)
                    * factorial(angular - k)
                    * factorial(angular - m - 2 * k)
                )
            )
            for y in range(m + 1):
                for a in range(k + 1):
                    for b in range(k - a + 1):
                        c = k - a - b
                        power = (m - y + 2 * a, y + 2 * b, angular - m - 2 * k + 2 * c)
                        multinomial = factorial(k) / (factorial(a) * factorial(b) * factorial(c))
                        polynomial[lookup[power]] += (
                            scale * legendre * comb(m, y) * 1j**y * multinomial
                        )
        matrix[angular + m] = polynomial.real
        if m:
            matrix[angular - m] = polynomial.imag
    powers.setflags(write=False)
    matrix.setflags(write=False)
    return powers, matrix


def _moments(alpha, beta, a, b, left_degree, right_degree):
    p = alpha + beta
    shift_a, shift_b = beta / p * (b - a), alpha / p * (a - b)
    table = np.zeros((left_degree + 1, right_degree + 1))
    table[0, 0] = np.sqrt(pi / p)
    for i in range(1, left_degree + 1):
        table[i, 0] = shift_a * table[i - 1, 0]
        if i > 1:
            table[i, 0] += (i - 1) / (2 * p) * table[i - 2, 0]
    for j in range(1, right_degree + 1):
        for i in range(left_degree + 1):
            table[i, j] = shift_b * table[i, j - 1]
            if i:
                table[i, j] += i / (2 * p) * table[i - 1, j - 1]
            if j > 1:
                table[i, j] += (j - 1) / (2 * p) * table[i, j - 2]
    return table


def radial_overlap(left: SphericalGaussian, right: SphericalGaussian, radius=None, *, order=None):
    """Equal-l radial block, with radius measured from the shared shell center."""
    _shell(left)
    _shell(right)
    if left.angular != right.angular or np.any(np.asarray(left.center) != right.center):
        raise ValueError("Same center and angular degree required")
    if radius is not None and (not np.isfinite(radius) or radius <= 0):
        raise ValueError("Positive finite radius required")
    if order is not None:
        if radius is None or type(order) is not int or not 2 <= order <= 512:
            raise ValueError("Positive sphere and bounded radial quadrature order required")
        nodes, weights = leggauss(order)
        r, weights = (nodes + 1) * radius / 2, weights * radius / 2
        values = [
            sum(c * np.exp(-e * r**2) for e, c in zip(s.exponents, s.amplitudes, strict=True))
            for s in (left, right)
        ]
        return (
            4
            * pi
            / (2 * left.angular + 1)
            * np.sum(weights * r ** (2 * left.angular + 2) * values[0].conj() * values[1])
        )
    h = left.angular + 1.5
    value = 0j
    for alpha, ca in zip(left.exponents, left.amplitudes, strict=True):
        for beta, cb in zip(right.exponents, right.amplitudes, strict=True):
            p = alpha + beta
            fraction = 1 if radius is None else gammainc(h, p * radius**2)
            value += np.conj(ca) * cb * 2 * pi / (2 * left.angular + 1) * gamma(h) / p**h * fraction
    return value


def shell_overlap(left: SphericalGaussian, right: SphericalGaussian, *, method="moments", extra=0):
    """All-space shell cross Gram by analytic moments or independent Gauss-Hermite."""
    _shell(left)
    _shell(right)
    if method not in ("moments", "hermite") or type(extra) is not int or not 0 <= extra <= 2:
        raise ValueError("Known method and at most two extra Hermite orders required")
    la, lb = left.angular, right.angular
    work = ((la + lb + 2) // 2 + extra) ** 3 * (2 * la + 2 * lb + 2)
    if work * len(left.exponents) * len(right.exponents) > 100_000_000:
        raise ValueError("Shell-pair point-function work cap exceeded")
    result = np.zeros((2 * la + 1, 2 * lb + 1), complex)
    a, b = np.array(left.center), np.array(right.center)
    if np.array_equal(a, b) and method == "moments":
        if la == lb:
            result = np.eye(2 * la + 1) * radial_overlap(left, right)
        return result
    pa, ma = _polynomials(la)
    pb, mb = _polynomials(lb)
    if method == "hermite":
        order = (la + lb + 2) // 2 + extra
        nodes, weights = hermgauss(order)
        coordinates = np.stack(np.meshgrid(nodes, nodes, nodes, indexing="ij"), axis=-1).reshape(
            -1, 3
        )
        w = np.prod(
            np.stack(np.meshgrid(weights, weights, weights, indexing="ij"), axis=-1), axis=-1
        ).ravel()
    for alpha, ca in zip(left.exponents, left.amplitudes, strict=True):
        for beta, cb in zip(right.exponents, right.amplitudes, strict=True):
            p = alpha + beta
            prefactor = np.conj(ca) * cb * np.exp(-alpha * beta / p * np.sum((a - b) ** 2))
            if method == "moments":
                cartesian = np.ones((len(pa), len(pb)))
                for d in range(3):
                    table = _moments(alpha, beta, a[d], b[d], la, lb)
                    cartesian *= table[pa[:, d, None], pb[None, :, d]]
                block = ma @ cartesian @ mb.T
            else:
                points = coordinates / np.sqrt(p) + (alpha * a + beta * b) / p
                va, vb = solid_harmonics(la, points - a), solid_harmonics(lb, points - b)
                block = va.T @ (w[:, None] * vb) / p**1.5
            result += prefactor * block
    if not np.isfinite(result).all():
        raise ValueError("Nonfinite Gaussian product integral")
    return result


def all_space_metric(shells, *, method="moments", extra=0, max_work=100_000_000):
    """Bounded all-space AO self Gram; Hermite work checked before matrix allocation."""
    shells, offsets = _basis(shells)
    if method not in ("moments", "hermite") or type(extra) is not int or not 0 <= extra <= 2:
        raise ValueError("Known method and at most two extra Hermite orders required")
    if type(max_work) is not int or not 1 <= max_work <= 100_000_000:
        raise ValueError("Positive bounded work cap required")
    work = sum(
        ((a.angular + b.angular + 2) // 2 + extra) ** 3
        * (2 * a.angular + 2 * b.angular + 2)
        * len(a.exponents)
        * len(b.exponents)
        for i, a in enumerate(shells)
        for b in shells[i:]
    )
    if work > max_work:
        raise ValueError("Gaussian-product point-function work cap exceeded")
    result = np.empty((offsets[-1], offsets[-1]), complex)
    for i, left in enumerate(shells):
        rows = slice(offsets[i], offsets[i + 1])
        for j in range(i, len(shells)):
            columns = slice(offsets[j], offsets[j + 1])
            block = shell_overlap(left, shells[j], method=method, extra=extra)
            result[rows, columns] = block
            if i != j:
                result[columns, rows] = block.conj().T
    return result


def generator_scales(shells, radius):
    """Tagged AO-row scales restoring origin-centered continuum norms at fixed a0."""
    shells, _ = _basis(shells)
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("Positive generation radius required")
    scales = []
    for shell in shells:
        value = radial_overlap(shell, shell, radius) if shell.continuum else 1.0
        if abs(np.imag(value)) > 1e-10 * max(abs(value), 1e-30):
            raise ValueError("Sphere self norm has an unexplained imaginary part")
        norm = float(np.real(value))
        if not np.isfinite(norm) or norm <= 0 or norm < 1e-12:
            raise ValueError("Continuum sphere norm is nonpositive or too ill-conditioned")
        scales.extend([1 / np.sqrt(norm)] * (2 * shell.angular + 1))
    return np.array(scales)


def uniform_ball_witness(shells, all_space, radius, scales):
    """Return restored hybrid Gram and conservative bound on its uniform-ball error.

    Only continuum-continuum entries are replaced by exact ball integrals.
    All other exterior integrals are bounded by contraction-aware Cauchy bounds.
    This is not zeroing unmeasured tails or a claim about an outgoing state's norm.
    """
    shells, offsets = _basis(shells)
    matrix, scales = np.asarray(all_space, complex), np.asarray(scales)
    n = offsets[-1]
    if (
        matrix.shape != (n, n)
        or scales.shape != (n,)
        or not np.isfinite(matrix).all()
        or not np.isrealobj(scales)
        or not np.isfinite(scales).all()
        or np.any(scales <= 0)
        or not np.isfinite(radius)
        or radius <= 0
    ):
        raise ValueError("Matching finite Gram, positive AO scales and sphere required")
    hybrid = matrix.copy()
    tail_roots, continuum = [], []
    for i, shell in enumerate(shells):
        h, degree = shell.angular + 1.5, 2 * shell.angular + 1
        lower = max(0.0, radius - np.linalg.norm(shell.center))
        root = sum(
            abs(c)
            * np.sqrt(2 * pi / degree * gamma(h) / (2 * e) ** h * gammaincc(h, 2 * e * lower**2))
            for e, c in zip(shell.exponents, shell.amplitudes, strict=True)
        )
        tail_roots.extend([root] * degree)
        continuum.extend([shell.continuum] * degree)
        if not shell.continuum:
            continue
        for j in range(i, len(shells)):
            other = shells[j]
            if not other.continuum:
                continue
            block = np.zeros((degree, 2 * other.angular + 1), complex)
            if shell.angular == other.angular:
                block = np.eye(degree) * radial_overlap(shell, other, radius)
            rows, columns = slice(offsets[i], offsets[i + 1]), slice(offsets[j], offsets[j + 1])
            hybrid[rows, columns], hybrid[columns, rows] = block, block.conj().T
    tail_roots = np.array(tail_roots) * scales
    error = np.outer(tail_roots, tail_roots)
    mask = np.logical_and.outer(continuum, continuum)
    error[mask] = 0.0
    return scales[:, None] * hybrid * scales[None, :], error


def shells_from_export(export):
    """Copy post-read shell amplitudes without interpreting them as raw file norms."""
    shells = tuple(
        SphericalGaussian(
            tuple(s["center"]),
            int(s["angular"]),
            tuple(s["exponents"]),
            tuple(s["normalization"] * s["contractions"] * s["primitive_norms"]),
            bool(s["continuum"]),
        )
        for s in export["shells"]
    )
    _, offsets = _basis(shells)
    descriptors = np.asarray(export["shell_descriptors"])
    if descriptors.shape != (6, len(shells)):
        raise ValueError("Complete CGTO-only shell descriptors required")
    for index in range(len(shells)):
        columns = descriptors[:, descriptors[1] == index + 1]
        if columns.shape != (6, 1) or columns[0, 0] != 1 or columns[3, 0] != offsets[index] + 1:
            raise ValueError("Unsupported nonsequential native AO/shell addressing")
    return shells
