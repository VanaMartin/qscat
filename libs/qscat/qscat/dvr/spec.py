"""Grid specification dataclasses for the FEM-DVR-ECS radial grid.

See `docs/physics/femdvr-ecs.md` for the construction this mirrors (ported from
eMoScat's FemDvrEcsGrid.cpp; the method is Rescigno & McCurdy, Phys. Rev. A 62,
032706 (2000)).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from qscat.exceptions import GridError

__all__ = ["ElementSpec", "GridSpec"]


@dataclass
class ElementSpec:
    """One finite element of the radial grid.

    `length` is the (real, positive) physical length of the element along the
    unscaled radial coordinate. `angle_deg` is the exterior-complex-scaling
    rotation angle (degrees); 0.0 means a real (unscaled) element.
    """

    length: float
    angle_deg: float = 0.0


@dataclass
class GridSpec:
    """Full grid specification: quadrature order (shared by all elements),
    the ordered list of elements, and the inner boundary x_min.

    Validates that complex (ECS) elements form a contiguous tail at the end
    of the element list -- eMoScat's `complex_negative` (ECS at the inner
    boundary) is deliberately not supported here.
    Computes the ECS pivot `R0 = x_min + sum(real element lengths)`, which by
    construction sits exactly on an element boundary.

    Caveat: using multiple *different* nonzero `angle_deg` values across tail
    elements (a bent/graded ECS contour) is REJECTED because it has never
    been validated here -- the validated, actually-used case is a single ECS
    tail angle shared by all complex elements.
    """

    quadrature: int
    elements: list[ElementSpec]
    x_min: float = 0.0
    R0: float = field(init=False)

    def __post_init__(self) -> None:
        """Validate the spec (quadrature >= 2, ordered elements) and derive `R0`."""
        if self.quadrature < 2:
            raise GridError("quadrature must be >= 2")
        if not self.elements:
            raise GridError("elements must be non-empty")

        seen_complex = False
        real_length_sum = 0.0
        distinct_nonzero_angles: set[float] = set()
        for el in self.elements:
            if el.length <= 0.0:
                raise GridError("element length must be positive")
            if el.angle_deg != 0.0:
                seen_complex = True
                distinct_nonzero_angles.add(el.angle_deg)
            else:
                if seen_complex:
                    raise GridError(
                        "complex (ECS) elements must be contiguous at the end of "
                        "the element list; found a real element after a complex one"
                    )
                real_length_sum += el.length

        if len(distinct_nonzero_angles) > 1:
            raise GridError(
                "a bent/graded ECS tail (more than one distinct nonzero "
                "angle_deg among the elements) is not validated and is "
                "rejected; use a single uniform ECS tail angle"
            )

        self.R0 = self.x_min + real_length_sum
