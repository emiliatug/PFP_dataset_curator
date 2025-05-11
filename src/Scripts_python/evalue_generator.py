from Modules.initial_preprocessing import InitialPreprocessor
from Utils.utils import *


def mainprocess(
    four_dict_loc: List[str],
    dict_dir_to_write_too: List[str],
    output_log_vars: List[str],
) -> None:
    output_dir = output_log_vars[0]
    os.makedirs(output_dir, exist_ok=True)

    # Part0.5 Set up logger file
    logger = generate_logging_file(
        name=f"{output_log_vars[1]}", loc_to_write=output_dir
    )
    logger.info("Starting the Global Part1: Data Generation: Part2")

    preprocessor = InitialPreprocessor()

    # Step1: Reading the 4 dict files
    df_main_all_go = load_parquet_files(four_dict_loc[0])
    df_main_all_ass = load_parquet_files(four_dict_loc[1])
    df_main_all_cp_go = load_parquet_files(four_dict_loc[2])
    df_main_all_cp_ass = load_parquet_files(four_dict_loc[3])
    df_child_parent = load_parquet_files(four_dict_loc[4])
    logger.info("Step 1 is done")

    # Step2: Reorganizing the dict of df_main_all_cp_go, df_main_all_cp_ass to make another 2 classes
    df_main_all_cp_go_perc = preprocessor.fast_merge_per(
        df_main_all_cp_go, df_child_parent
    )
    df_main_all_cp_ass_perc = preprocessor.fast_merge_per(
        df_main_all_cp_ass, df_child_parent
    )
    logger.info("Step 2 is done")

    # Step3: filters the 4 dict (that are now df, by removing the rows where the query protein is the same as the protein where it was found)
    df_main_all_go_filter, _ = preprocessor.filter_by_self(
        df_main_all_go, ["Out_of_Total", "Verbal_Desc", "Type_of_GO"]
    )
    df_main_all_ass_filter, _ = preprocessor.filter_by_self(
        df_main_all_ass, ["Out_of_Total", "Verbal_Desc", "Type_of_GO"]
    )
    df_main_all_cp_go_perc_filter, _ = preprocessor.filter_by_self(
        df_main_all_cp_go_perc,
        ["Out_of_Total", "Verbal_Desc", "GO_Term_Child", "Type_of_GO", "Parent"],
    )
    df_main_all_cp_ass_perc_filter, _ = preprocessor.filter_by_self(
        df_main_all_cp_ass_perc,
        [
            "Out_of_Total",
            "Verbal_Desc",
            "GO_Term_Child",
            "Type_of_GO",
            "found_or_ass_go",
            "Parent",
        ],
    )
    logger.info("Step 3 is done")

    # Step4: calculates the -log(Evalue)+b
    df_main_all_go_hist0 = preprocessor.apply_neglog10_transform(df_main_all_go_filter)
    df_main_all_ass_hist0 = preprocessor.apply_neglog10_transform(
        df_main_all_ass_filter
    )
    df_main_all_cp_go_hist0 = preprocessor.apply_neglog10_transform(
        df_main_all_cp_go_perc_filter
    )
    df_main_all_cp_ass_hist0 = preprocessor.apply_neglog10_transform(
        df_main_all_cp_ass_perc_filter
    )
    logger.info("Step 4 is done")

    # Step5: filters the 5: saves the files
    saving_lists_of_locats_and_dicts(
        dict_dir_to_write_too,
        [
            df_main_all_go_hist0,
            df_main_all_ass_hist0,
            df_main_all_cp_go_hist0,
            df_main_all_cp_ass_hist0,
        ],
        logger,
    )
    logger.info("hist file is saved")

    logger.info("Step 5 is done")
    logger.info("Global1_part2.py is done")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python script.py <four_dict_loc>, <dict_dir_to_write_too>,<output_log_vars>"
        )
        sys.exit(1)

    four_dict_loc = [col.strip() for col in sys.argv[1].split(",")]
    dict_dir_to_write_too = [col.strip() for col in sys.argv[2].split(",")]
    output_log_vars = [col.strip() for col in sys.argv[3].split(",")]
    mainprocess(four_dict_loc, dict_dir_to_write_too, output_log_vars)
