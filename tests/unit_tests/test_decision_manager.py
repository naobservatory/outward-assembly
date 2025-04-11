import pytest
import yaml
from outward_assembly.decision_manager import DecisionManager
from outward_assembly.strategy_helper import (
    CONTIG_COUNT,
    LONGEST_CONTIG,
    READ_PAIR_COUNT,
)
from outward_assembly.actions import decrease_k, next_priority, increase_k

# Sample YAML configuration for testing
YAML_CONFIG = """
decision:
  automate: true
  strategy: example_strategy
"""


@pytest.mark.parametrize(
    "outward_assembly_metrics, expected_actions",
    [
        (
            {READ_PAIR_COUNT: 90, LONGEST_CONTIG: 150, CONTIG_COUNT: 3},
            [lambda metrics: next_priority(metrics)],
        ),
        (
            {READ_PAIR_COUNT: 90, LONGEST_CONTIG: 150, CONTIG_COUNT: 1},
            [
                lambda metrics: next_priority(metrics),
                lambda metrics: decrease_k(3)(metrics),
            ],
        ),
        (
            {READ_PAIR_COUNT: 100, LONGEST_CONTIG: 200, CONTIG_COUNT: 1},
            [lambda metrics: decrease_k(3)(metrics)],
        ),
        (
            {READ_PAIR_COUNT: 1001, LONGEST_CONTIG: 200, CONTIG_COUNT: 2},
            [lambda metrics: increase_k(5)(metrics)],
        ),
    ],
)
def test_decision_manager_actions(outward_assembly_metrics, expected_actions):
    """Test the DecisionManager's ability to evaluate metrics and return correct actions."""
    config = yaml.safe_load(YAML_CONFIG)
    manager = DecisionManager(config["decision"])
    actions = manager.collect_actions(outward_assembly_metrics)

    # Check that we got the expected number of actions
    assert len(actions) == len(expected_actions)
