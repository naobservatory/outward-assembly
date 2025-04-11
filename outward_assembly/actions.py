from typing import Dict, Any, Callable


def increase_k(val: int) -> Callable[[Dict[str, Any]], int]:
    """Creates an action that increases k by the specified value.
    Returns a function that computes the new k value given a state."""
    return lambda state: dict(state, read_subset_k=state["read_subset_k"] + val)


def decrease_k(val: int) -> Callable[[Dict[str, Any]], int]:
    """Creates an action that decreases k by the specified value.
    Returns a function that computes the new k value given a state."""
    return lambda state: dict(state, read_subset_k=state["read_subset_k"] - val)


def next_priority() -> Callable[[Dict[str, Any]], int]:
    """Creates an action that advances to the next priority.
    Returns a function that computes the new priority value given a state."""
    return lambda state: dict(state, dataset_priority=state["dataset_priority"] + 1)
