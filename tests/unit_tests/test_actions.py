from outward_assembly.actions import increase_k, decrease_k, next_priority


def test_increase_k():
    """Test the increase_k action."""
    curr_exec_state = {"read_subset_k": 10}
    curr_exec_state = increase_k(5)(
        curr_exec_state
    )  # Call the returned function with the state
    assert curr_exec_state["read_subset_k"] == 15


def test_decrease_k():
    """Test the decrease_k action."""
    curr_exec_state = {"read_subset_k": 10}
    curr_exec_state = decrease_k(3)(
        curr_exec_state
    )  # Call the returned function with the state
    assert curr_exec_state["read_subset_k"] == 7


def test_next_priority():
    """Test the next_priority action."""
    curr_exec_state = {"dataset_priority": 1}
    curr_exec_state = next_priority()(
        curr_exec_state
    )  # Call the returned function with the state
    assert curr_exec_state["dataset_priority"] == 2
