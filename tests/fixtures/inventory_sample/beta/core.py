"""Beta core module."""

from ..alpha.core import shared_helper


def documented_public(a, b):
    """A homonym of alpha.core.documented_public with a different body."""
    return a - b


def clone_b(items):
    """Count truthy items."""
    n = 0
    for item in items:
        if item:
            n += 1
    return n


def near_clone_b(items, limit):
    """Sum the items at or below limit, and count the ones skipped."""
    total = 0
    skipped = 0
    for item in items:
        if item <= limit:
            total += item
        else:
            skipped += 1
    return total, skipped


def uses_helper(values):
    """Call across packages."""
    return shared_helper(values)
