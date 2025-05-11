from Modules.train_val_generator_handler import TrainValidationSplitHandler
from Utils.utils import *


def mainprocess(
    input1: List[str],
    output_dir_files: str,
    output_log_vars: List[str],
) -> None:
    output_dir = output_log_vars[0]
    os.makedirs(output_dir, exist_ok=True)

    logger = generate_logging_file(
        name=f"{output_log_vars[1]}", loc_to_write=output_dir
    )
    logger.info(
        "Starting the Global Part1: Data Generation: Part3-Splitting entire raw data into 5 CV"
    )

    """This fx reads the 4 dict, makes them a single dictionary, and split them into test and train and val CV; saving them at the end."""

    split_handler = TrainValidationSplitHandler()

    all_dictionaries = split_handler.load_4_class_dictionaries(
        input1[0], input1[1], input1[2], input1[3]
    )

    # generated test/train/validations
    combined_df = split_handler.combine_4_class_dictionaries(
        all_dictionaries["df_main_allGO_hist0"],
        all_dictionaries["df_main_allAss_hist0"],
        all_dictionaries["df_main_allCP_GO_hist0"],
        all_dictionaries["df_main_allCP_Ass_hist0"],
    )

    logger.info("generating_combined_df_from_4_main_dict is done")

    split_handler.perform_group_shuffle_split(
        logger,
        combined_df,
        output_dir_files,
        n_splits_hold_out=1,
        n_splits_outer=5,
        n_splits_inner=1,
    )

    logger.info("generating_test_train_and_5CV_and_processing_them is done")
    logger.info("Process is complete")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input1> <input2> <input3>,<output_log_vars>")
        sys.exit(1)

    output_dir_files = sys.argv[2]  # Collect all arguments after the script name
    input1 = [col.strip() for col in sys.argv[1].split(",")]
    output_log_vars = [col.strip() for col in sys.argv[3].split(",")]

    mainprocess(input1, output_dir_files, output_log_vars)
