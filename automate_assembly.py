from outward_assembly import outward_assembly
from outward_assembly.io_helpers import (
    get_s3_paths_by_priority,
    load_config,
)
from outward_assembly.execution_state import ExecutionState
from outward_assembly.decision_manager import DecisionManager
from dotenv import load_dotenv
import logging
import sys
import os

load_dotenv()


def main(input_config_path: str):
    """
    Main function for setting up the classes necessary for running outward assembly in a loop based on the user defined parameters.

    See docs/usage.md for more information.
    """
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    config = load_config(input_config_path)

    execution_state = ExecutionState.from_config(config)

    manager = DecisionManager(config["decision"])

    final_output_path = (
        execution_state.output_dir + "/" + execution_state.output_filename
    )

    while True:
        input_data = get_s3_paths_by_priority(
            execution_state.input_dataset_list, execution_state.dataset_priority
        )

        if input_data == []:
            logger.info("Pipeline complete: no more datasets to assemble")
            break

        output_path = (
            execution_state.output_dir
            + "/"
            + "outer_iter_"
            + str(execution_state.current_outer_iterations)
            + "_"
            + execution_state.output_filename
        )

        assembly_metrics = outward_assembly(
            s3_paths=input_data,
            seed_path=execution_state.input_seed_path,
            output_path=output_path,
            adapters_path=execution_state.adapter_path,
            work_dir_parent=execution_state.work_dir,
            read_subset_k=execution_state.read_subset_k,
            use_batch=execution_state.use_batch,
        )

        # @TODO: Implement contig analysis as another module, and then add those metrics everywhere relevant

        actions_to_run = manager.collect_actions(assembly_metrics)

        current_state = execution_state.get_adjustable_parameters()
        for action in actions_to_run:
            current_state = action(current_state)

        execution_state.update_adjustable_parameters(current_state)

        # Write the state to a file so that it can be resumed from the same point after any changes were specified
        work_dir = assembly_metrics["work_dir"]
        execution_state.write_config(work_dir + "/config.yaml")

        if not manager.automate:
            logger.info("Pipeline complete: automation disabled")
            os.symlink(output_path, final_output_path)
            break

        if execution_state.check_limits():
            logger.info("Pipeline complete: limits reached")
            os.symlink(output_path, final_output_path)
            break

        if actions_to_run == []:
            logger.info("Pipeline complete: actions are not met")
            os.symlink(output_path, final_output_path)
            break


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: The input configuration path is required.")
        print("Variable: input_config - Path to the input configuration YAML file.")
        sys.exit(1)

    input_config_path = sys.argv[1]
    main(input_config_path)
