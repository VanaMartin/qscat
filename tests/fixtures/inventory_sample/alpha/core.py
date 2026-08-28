"""Alpha core module."""


def documented_public(a, b):
    """Add two numbers and return the sum."""
    return a + b


def shared_helper(values):
    """Return the mean of values."""
    total = 0
    for v in values:
        total += v
    return total / len(values)


def clone_a(items):
    """Count truthy items."""
    n = 0
    for item in items:
        if item:
            n += 1
    return n


def near_clone_a(items, limit):
    """Sum the items below limit, and count the ones skipped."""
    total = 0
    skipped = 0
    for item in items:
        if item < limit:
            total += item
        else:
            skipped += 1
    return total, skipped


def _private_dead(x):
    """Nothing calls this."""
    return x * 2


def string_dispatched(x):
    """Reachable only through config/sample.yaml."""
    return x + 1


def wide(a, b, c, d, e, f, g, h, i):
    """Nine parameters, no branching."""
    return (a, b, c, d, e, f, g, h, i)


def undocumented_export(x):
    return x
