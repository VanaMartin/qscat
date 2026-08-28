"""One comment of each rubric category."""


def stepper(dt, order):
    """Return a propagation coefficient."""
    # loop over the orders
    total = 0.0
    for k in range(order):
        # we originally used a dense solve here, but it was too slow
        # PRA 47 Eq. (2.15): the local limit drops the memory integral
        total += dt**k
    # do not reorder: the caller relies on the k=0 term arriving first
    return total
