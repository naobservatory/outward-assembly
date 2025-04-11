"""
Helper module that defines the variables on which the user can condition their strategy.
"""

from typing import TypedDict, Literal

# Define the metric keys as string literals
CONTIG_COUNT: Literal["final_contig_count"] = "final_contig_count"
LONGEST_CONTIG: Literal["final_longest_contig_length"] = "final_longest_contig_length"
READ_PAIR_COUNT: Literal["final_read_pair_counts"] = "final_read_pair_counts"


class OutwardAssemblyMetrics(TypedDict):
    """TypedDict defining the structure of metrics passed to strategy functions"""

    final_contig_count: int
    final_longest_contig_length: int
    final_read_pair_counts: int
