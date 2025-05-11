import gc

import matplotlib.pyplot as plt
import numpy as np

from Modules.final_dataset_handler import FinalDatasetGenerationHandler
from Utils.utils import *


def main_process(
    prediction_inputs: List[str],
    sparse_inputs: List[str],
    run_metadata: List[str],
) -> None:
    """
    Main process for Global Part1: Part6.
    Generates the final training dataset 'final_train_df.pkl' from prediction and sparse array inputs.

    Parameters:
    prediction_inputs (List[str]): Paths to prediction input files.
    sparse_inputs (List[str]): Paths to sparse arrays and bin edge inputs.
    run_metadata (List[str]): Run-specific metadata including output directories and run identifiers.
    """
    output_dir = run_metadata[2]
    os.makedirs(output_dir, exist_ok=True)

    logger = generate_logging_file(
        name=f"{run_metadata[3]}_{run_metadata[1]}", loc_to_write=output_dir
    )
    logger.info(
        "Starting the Global Part1: Part6: Generating the Final Train/Train_val.final_train_df.pkl"
    )

    if not all(os.path.exists(path) for path in prediction_inputs):
        raise ValueError(
            "One or more input directories in `prediction_inputs` do not exist."
        )

    single_prediction_with_depth = pd.read_pickle(prediction_inputs[0])
    sparse_aa = pd.read_pickle(sparse_inputs[0])
    sparse_bb = pd.read_pickle(sparse_inputs[1])
    sparse_cc = pd.read_pickle(sparse_inputs[2])
    sparse_dd = pd.read_pickle(sparse_inputs[3])

    bin_edges_aa = np.load(sparse_inputs[4])
    bin_edges_bb = np.load(sparse_inputs[5])
    bin_edges_cc = np.load(sparse_inputs[6])
    bin_edges_dd = np.load(sparse_inputs[7])

    output_dir_files = f"{run_metadata[0]}/{run_metadata[1]}"
    os.makedirs(output_dir_files, exist_ok=True)

    scores_df = single_prediction_with_depth[
        [
            "Protein",
            "Raw_Score",
            "Original_PFP_Pred_GO_term",
            "BinaryScore",
            "LinScore",
            "Root_Flag",
            "Depth",
        ]
    ]

    handler = FinalDatasetGenerationHandler()

    sparse_aa = handler.read_clean_part4_var(sparse_aa)
    sparse_bb = handler.read_clean_part4_var(sparse_bb)
    sparse_cc = handler.read_clean_part4_var(sparse_cc)
    sparse_dd = handler.read_clean_part4_var(sparse_dd)
    print("Step1 is done")
    logger.info("Step1 is done")

    scores_counts_aa = handler.combine_scores_and_counts(scores_df, sparse_aa)
    scores_counts_bb = handler.combine_scores_and_counts(scores_df, sparse_bb)
    scores_counts_cc = handler.combine_scores_and_counts(scores_df, sparse_cc)
    scores_counts_dd = handler.combine_scores_and_counts(scores_df, sparse_dd)
    print("Step2 is done")
    logger.info("Step2 is done")

    midpoints_aa = handler.find_midpoints(list(sparse_aa.columns))
    midpoints_bb = handler.find_midpoints(list(sparse_bb.columns))
    midpoints_cc = handler.find_midpoints(list(sparse_cc.columns))
    midpoints_dd = handler.find_midpoints(list(sparse_dd.columns))
    print("Step3 is done")
    logger.info("Step3 is done")

    new_ranges_aa, range_list_aa = handler.make_ranges_from_bins(
        midpoints_aa, bin_edges_aa
    )
    new_ranges_bb, range_list_bb = handler.make_ranges_from_bins(
        midpoints_bb, bin_edges_bb
    )
    new_ranges_cc, range_list_cc = handler.make_ranges_from_bins(
        midpoints_cc, bin_edges_cc
    )
    new_ranges_dd, range_list_dd = handler.make_ranges_from_bins(
        midpoints_dd, bin_edges_dd
    )

    range_list_aa1 = [x + ["AA"] for x in range_list_aa]
    range_list_bb1 = [x + ["BB"] for x in range_list_bb]
    range_list_cc1 = [x + ["CC"] for x in range_list_cc]
    range_list_dd1 = [x + ["DD"] for x in range_list_dd]
    print("Step4 is done")
    logger.info("Step4 is done")

    renamed_aa = handler.rename_df(scores_counts_aa, range_list_aa1)
    renamed_bb = handler.rename_df(scores_counts_bb, range_list_bb1)
    renamed_cc = handler.rename_df(scores_counts_cc, range_list_cc1)
    renamed_dd = handler.rename_df(scores_counts_dd, range_list_dd1)
    print("Step5 is done")
    logger.info("Step5 is done")

    tmp1 = pd.merge(
        renamed_aa,
        renamed_bb,
        on=["Name", "Raw_Score", "GO", "LinScore", "Root_Flag", "Depth"],
        how="outer",
    )
    tmp2 = pd.merge(
        tmp1,
        renamed_cc,
        on=["Name", "Raw_Score", "GO", "LinScore", "Root_Flag", "Depth"],
        how="outer",
    )
    tmp3 = pd.merge(
        tmp2,
        renamed_dd,
        on=["Name", "Raw_Score", "GO", "LinScore", "Root_Flag", "Depth"],
        how="outer",
    )
    tmp3 = tmp3.fillna(0)

    tmp4 = tmp3.sort_values("Raw_Score").reset_index(drop=True)
    tmp4_to_analyze_df = tmp4.copy()
    tmp4_to_analyze_df.to_pickle(
        f"./{output_dir_files}/tmp4_lin_to_analyze_df.pkl", protocol=4
    )
    print("Step6 is done")
    logger.info("Step6 is done")

    bgraph_aa = handler.make_ranges_for_bgraph(range_list_aa)
    bgraph_dd = handler.make_ranges_for_bgraph(range_list_dd)
    print("Step7 is done")
    logger.info("Step7 is done")

    colors = plt.get_cmap("hsv", 100)
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 4))
    handler.plot_ranges(ax1, bgraph_aa, colors, "Bins from Class A")
    handler.plot_ranges(ax2, bgraph_dd, colors, "Bins from Class D")
    print("Step8 is done")
    logger.info("Step8 is done")

    plt.tight_layout()
    fig.savefig(f"{output_dir_files}/my_plot.png", dpi=100)

    single_prediction_with_depth = single_prediction_with_depth.drop(
        [
            "DirectBlastAA",
            "AssociatedBlastBB",
            "DirectParentCC",
            "AssociatedParentDD",
            "PFP_Related/NearestNeighbor",
            "Root_Flag",
        ],
        axis=1,
    )

    if single_prediction_with_depth["Protein"].equals(
        single_prediction_with_depth["Name"]
    ):
        single_prediction_with_depth = single_prediction_with_depth.drop(
            "Protein", axis=1
        )

    tmp4_to_analyze_df.rename(columns={"GO": "Original_PFP_Pred_GO_term"}, inplace=True)
    merged_df = pd.merge(
        single_prediction_with_depth,
        tmp4_to_analyze_df,
        on=["Name", "Raw_Score", "Original_PFP_Pred_GO_term", "LinScore", "Depth"],
        how="inner",
    )
    merged_df1 = merged_df.drop(
        ["Predicted_Ontology_Retained", "Predicted_Ontology_C", "Root_Flag"], axis=1
    )
    merged_df1 = merged_df1.drop_duplicates()

    col_l = [
        "Feat_Freq" + str(x) for x in range(0, len(merged_df1.columns.to_list()[14:]))
    ]
    merged_df1.columns.values[14:] = col_l
    merged_df1.columns = merged_df1.columns.str.strip()
    merged_df.to_pickle(
        f"./{output_dir_files}/pre_named_final_train_df.pkl", protocol=4
    )
    merged_df1.to_pickle(f"./{output_dir_files}/final_train_df.pkl", protocol=4)
    logger.info("Step9 is done")
    print("Finished the main script", flush=True)
    logger.info("Global1_part6.py is finished")
    gc.collect()
    logging.shutdown()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python script.py <prediction_inputs> <sparse_inputs> <run_metadata>"
        )
        sys.exit(1)

    prediction_inputs = [col.strip() for col in sys.argv[1].split(",")]
    sparse_inputs = [col.strip() for col in sys.argv[2].split(",")]
    run_metadata = [col.strip() for col in sys.argv[3].split(",")]
    main_process(prediction_inputs, sparse_inputs, run_metadata)
