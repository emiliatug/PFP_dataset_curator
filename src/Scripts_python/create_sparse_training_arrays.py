import gc

import numpy as np

from Modules.sparse_training_array_handler import SparseTrainingArrayHandler
from Modules.train_val_generator_handler import TrainValidationSplitHandler
from Utils.utils import *


def mainprocess(
    class_data_paths: List[str],
    training_data_and_predictions: List[str],
    output_and_run_metadata: List[str],
) -> None:
    output_dir = output_and_run_metadata[2]
    os.makedirs(output_dir, exist_ok=True)

    logger = generate_logging_file(
        name=f"{output_and_run_metadata[3]}_{output_and_run_metadata[1]}",
        loc_to_write=output_dir,
    )
    logger.info("Starting the Global Part1: Data Generation:Train: Part4")

    if not all(os.path.exists(path) for path in class_data_paths):
        raise ValueError(
            "One or more input directories in `class_data_paths` do not exist."
        )

    handler = SparseTrainingArrayHandler()

    cv_subset_training = pd.read_parquet(training_data_and_predictions[0])
    logger.info(f"This is the size of the training data {cv_subset_training.shape}")

    original_prediction_table_dir = training_data_and_predictions[1]
    output_run_dir = f"{output_and_run_metadata[0]}/{output_and_run_metadata[1]}"
    os.makedirs(output_run_dir, exist_ok=True)

    all_class_data_dicts = TrainValidationSplitHandler().load_4_class_dictionaries(
        class_data_paths[0],
        class_data_paths[1],
        class_data_paths[2],
        class_data_paths[3],
    )
    logger.info("Step1 is done")

    new_shape_before_hist_a = pd.concat(
        all_class_data_dicts["df_main_all_go_hist0"].values()
    )
    new_shape_before_hist_b = pd.concat(
        all_class_data_dicts["df_main_all_ass_hist0"].values()
    )
    new_shape_before_hist_c = pd.concat(
        all_class_data_dicts["df_main_all_cp_go_hist0"].values()
    )
    new_shape_before_hist_d = pd.concat(
        all_class_data_dicts["df_main_all_cp_ass_hist0"].values()
    )
    logger.info("Step2 is done")

    logger.info(
        f"Shape of the 4 classes from the original data is {new_shape_before_hist_a.shape}, {new_shape_before_hist_b.shape}, {new_shape_before_hist_c.shape}, {new_shape_before_hist_d.shape}"
    )

    class_a_cv_subset, class_b_cv_subset, class_c_cv_subset, class_d_cv_subset = (
        handler.split_dataframe_by_class(cv_subset_training)
    )

    class_a_enriched = handler.enrich_dataframe_by_class(
        class_a_cv_subset, new_shape_before_hist_a
    )
    class_b_enriched = handler.enrich_dataframe_by_class(
        class_b_cv_subset, new_shape_before_hist_b
    )
    class_c_enriched = handler.enrich_dataframe_by_class(
        class_c_cv_subset, new_shape_before_hist_c
    )
    class_d_enriched = handler.enrich_dataframe_by_class(
        class_d_cv_subset, new_shape_before_hist_d
    )
    logger.info("Step 2 & 3 is done")

    range_calculator = handler.RangeCalculator()
    range_a, range_e_a = range_calculator.give_range(class_a_enriched, "-log+b")
    range_b, range_e_b = range_calculator.give_range(class_b_enriched, "-log+b")
    range_c, range_e_c = range_calculator.give_range(class_c_enriched, "-log+b_P")
    range_d, range_e_d = range_calculator.give_range(class_d_enriched, "-log+b_P")
    range_calculator.reset_count()
    logger.info("Step 4 is done")

    df_a = handler.make_train_df_to_dic(class_a_enriched)
    df_b = handler.make_train_df_to_dic(class_b_enriched)
    df_c = handler.make_train_df_to_dic(class_c_enriched)
    df_d = handler.make_train_df_to_dic(class_d_enriched)
    logger.info("Step 5 is done")

    df_aa_sparse, edges_aa, midpoints_aa = handler.get_histogram_values_sparse(
        df_a, "-log+b", "auto", range_a, "_AA"
    )
    df_bb_sparse, edges_bb, midpoints_bb = handler.get_histogram_values_sparse(
        df_b, "-log+b", "auto", range_b, "_BB"
    )
    df_cc_sparse, edges_cc, midpoints_cc = handler.get_histogram_values_sparse(
        df_c, "-log+b_P", "auto", range_c, "_CC"
    )
    df_dd_sparse, edges_dd, midpoints_dd = handler.get_histogram_values_sparse(
        df_d, "-log+b_P", "auto", range_d, "_DD"
    )
    logger.info("Step 6 is done")

    df_aa_sparse.to_pickle(f"{output_run_dir}/df_aa_sparse_aray.pkl", protocol=4)
    df_bb_sparse.to_pickle(f"{output_run_dir}/df_bb_sparse_aray.pkl", protocol=4)
    df_cc_sparse.to_pickle(f"{output_run_dir}/df_cc_sparse_aray.pkl", protocol=4)
    df_dd_sparse.to_pickle(f"{output_run_dir}/df_dd_sparse_aray.pkl", protocol=4)

    np.save(f"{output_run_dir}/df_aa_bin_edges.npy", edges_aa)
    np.save(f"{output_run_dir}/df_bb_bin_edges.npy", edges_bb)
    np.save(f"{output_run_dir}/df_cc_bin_edges.npy", edges_cc)
    np.save(f"{output_run_dir}/df_dd_bin_edges.npy", edges_dd)

    np.save(f"{output_run_dir}/df_aa_bin_midpoints.npy", midpoints_aa)
    np.save(f"{output_run_dir}/df_bb_bin_midpoints.npy", midpoints_bb)
    np.save(f"{output_run_dir}/df_cc_bin_midpoints.npy", midpoints_cc)
    np.save(f"{output_run_dir}/df_dd_bin_midpoints.npy", midpoints_dd)
    logger.info("Step 7 is done")

    df_aa_filtered = handler.drop_duplicates_sparse(df_aa_sparse, "filtered_count_AA")
    df_bb_filtered = handler.drop_duplicates_sparse(df_bb_sparse, "filtered_count_BB")
    df_cc_filtered = handler.drop_duplicates_sparse(df_cc_sparse, "filtered_count_CC")
    df_dd_filtered = handler.drop_duplicates_sparse(df_dd_sparse, "filtered_count_DD")

    sparse_aa = handler.convert_sparse_column_to_dataframe(
        df_aa_filtered, "filtered_count_AA"
    )
    sparse_bb = handler.convert_sparse_column_to_dataframe(
        df_bb_filtered, "filtered_count_BB"
    )
    sparse_cc = handler.convert_sparse_column_to_dataframe(
        df_cc_filtered, "filtered_count_CC"
    )
    sparse_dd = handler.convert_sparse_column_to_dataframe(
        df_dd_filtered, "filtered_count_DD"
    )
    logger.info("Step 8 is done")

    class_a_merged_tmp = handler.merge_sparse_columns_fixed(sparse_aa, target_cols=100)
    class_a_merged_final = handler.add_sum_row(class_a_merged_tmp)
    class_b_merged_tmp = handler.merge_sparse_columns_fixed(sparse_bb, target_cols=50)
    class_b_merged_final = handler.add_sum_row(class_b_merged_tmp)
    class_d_merged_tmp = handler.merge_sparse_columns_fixed(sparse_dd, target_cols=100)
    class_d_merged_final = handler.add_sum_row(class_d_merged_tmp)

    class_c_merged_final = handler.group_sparse_matrix(
        sparse_cc.iloc[:, 2:], batch_size=400
    )
    logger.info("Step 9 & 10 is done")

    class_a_merged_final = class_a_merged_final.iloc[:-1, :]
    class_b_merged_final = class_b_merged_final.iloc[:-1, :]
    class_c_final = pd.concat(
        [sparse_cc.iloc[:, :2], class_c_merged_final.iloc[:-5, :]], axis=1
    )
    class_d_merged_final = class_d_merged_final.iloc[:-1, :]

    class_a_merged_final.to_pickle(f"{output_run_dir}/aa33.pkl", protocol=4)
    class_b_merged_final.to_pickle(f"{output_run_dir}/bb33.pkl", protocol=4)
    class_c_final.to_pickle(f"{output_run_dir}/cc33.pkl", protocol=4)
    class_d_merged_final.to_pickle(f"{output_run_dir}/dd33.pkl", protocol=4)
    logger.info("Step 11 is done")

    df_main_100 = handler.generate_original_predictions(original_prediction_table_dir)
    df_main_100_df = (
        pd.concat(df_main_100.values())
        .drop("Verbal_Desc", axis=1)
        .rename(columns={"GO_Term": "GO"})
    )

    tmp12 = pd.merge(
        df_aa_filtered, df_bb_filtered, on=["GO", "Name"], how="outer"
    ).drop("_merge", axis=1)
    tmp123 = pd.merge(tmp12, df_cc_filtered, on=["GO", "Name"], how="outer").drop(
        "_merge", axis=1
    )
    tmp1234 = pd.merge(tmp123, df_dd_filtered, on=["GO", "Name"], how="outer").drop(
        "_merge", axis=1
    )
    final_merged = pd.merge(df_main_100_df, tmp1234, on=["Name", "GO"], how="inner")

    final_merged.to_pickle(f"{output_run_dir}/jonedABCD_PFP.pkl", protocol=4)
    logger.info("Step 12 is done")

    logger.info("Global1_part4.py is finished")
    gc.collect()
    logging.shutdown()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python script.py <class_data_paths> <training_data_and_predictions> <output_and_run_metadata>"
        )
        sys.exit(1)

    class_data_paths = [col.strip() for col in sys.argv[1].split(",")]
    training_data_and_predictions = [col.strip() for col in sys.argv[2].split(",")]
    output_and_run_metadata = [col.strip() for col in sys.argv[3].split(",")]

    mainprocess(
        class_data_paths, training_data_and_predictions, output_and_run_metadata
    )
