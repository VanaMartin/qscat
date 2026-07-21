"""Grid specification dataclasses for the FEM-DVR-ECS radial grid.

See .superpowers/sdd/femdvr-ecs-extraction.md sections 1-3 and
docs/superpowers/specs/2026-07-21-femdvr-ecs-grid-design.md for the construction
this mirrors (ported from eMoScat's FemDvrEcsGrid.cpp).
"""

from dataclasses import dataclass, field


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
    boundary) is out of scope for this port (see design spec, Out of scope).
    Computes the ECS pivot `R0 = x_min + sum(real element lengths)`, which by
    construction sits exactly on an element boundary.

    Caveat: using multiple *different* nonzero `angle_deg` values across tail
    elements (a bent/graded ECS contour) is REJECTED UNTIL VALIDATED IN
    SUB-PROJECT #2 -- the validated, actually-used case is a single ECS tail
    angle shared by all complex elements.
    """

    quadrature: int
    elements: list[ElementSpec]
    x_min: float = 0.0
    R0: float = field(init=False)

    def __post_init__(self) -> None:
        if self.quadrature < 2:
            raise ValueError("quadrature must be >= 2")
        if not self.elements:
            raise ValueError("elements must be non-empty")

        seen_complex = False
        real_length_sum = 0.0
        distinct_nonzero_angles: set[float] = set()
        for el in self.elements:
            if el.length <= 0.0:
                raise ValueError("element length must be positive")
            if el.angle_deg != 0.0:
                seen_complex = True
                distinct_nonzero_angles.add(el.angle_deg)
            else:
                if seen_complex:
                    raise ValueError(
                        "complex (ECS) elements must be contiguous at the end of "
                        "the element list; found a real element after a complex one"
                    )
                real_length_sum += el.length

        if len(distinct_nonzero_angles) > 1:
            raise ValueError(
                "a bent/graded ECS tail (more than one distinct nonzero "
                "angle_deg among the elements) is rejected until validated "
                "in sub-project #2; use a single uniform ECS tail angle"
            )

        self.R0 = self.x_min + real_length_sum
