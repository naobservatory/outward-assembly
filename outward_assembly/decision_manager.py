# decision_manager.py
import logging
from typing import Any, Callable, Dict, List

from outward_assembly import strategy

logger = logging.getLogger(__name__)


class DecisionManager:
    """
    Import the user defined strategy, and manage the decisions that are made based on the strategy.
    """

    def __init__(self, config: Dict[str, Any]):
        self.automate = config.get("automate", False)
        self.strategy_name = config.get("strategy")
        # Check if the strategy exists in the strategy module's functions
        if self.strategy_name and getattr(strategy, self.strategy_name):
            self.strategy_fn = getattr(strategy, self.strategy_name)
        else:
            self.strategy_fn = None
            if self.automate:
                raise ValueError(
                    f"Cannot automate decisions without a valid 'strategy' in config. "
                    f"Found: {self.strategy_name}"
                )

    def collect_actions(
        self, outward_assembly_metrics: Dict[str, Any]
    ) -> List[Callable[[Dict[str, Any]], Any]]:
        """
        If the user has specified a strategy that exists and `automate` is True,
        return the list of actions to take.
        """
        if not self.automate or not self.strategy_fn:
            return []

        # Call the strategy function to get the actions to take
        return self.strategy_fn(outward_assembly_metrics)
