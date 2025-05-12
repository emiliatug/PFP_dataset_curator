import os
import os
import pickle
import sys
from collections import defaultdict
from functools import reduce

import pandas as pd

from Modules.pygosemsim_handler import PyGOsemsimHandler
from Modules.sparse_training_array_handler import SparseTrainingArrayHandler
from Modules.train_val_generator_handler import TrainValidationSplitHandler

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
import gc

from Modules.val_test_count_handler import *
from Utils.utils import *
from pygosemsim import (
    graph,
    similarity,
)

execution_dir = os.getcwd()
os.chdir(execution_dir)


def mainprocess(
    input1: List[str],
    input2: List[str],
    input3: List[str],
    org_var: List[str],
    input5: List[str],
    input6: List[str],
) -> None:

    output_dir = org_var[2]
    os.makedirs(output_dir, exist_ok=True)

    # Part0.5 Set up logger file
    logger = generate_logging_file(
        name=f"{org_var[3]}_{org_var[1]}", loc_to_write=output_dir
    )
    logger.info("Starting the Global Part1: Part6: either val or train")

    # Step0: Real most of the the inputs
    print("Most of the inputs are: ")
    print("Name1.0", "val_or_test_df", input1[0])

    print("Name1.1", "updated_uniprot_6_species_with_goterms_gaf", input1[1])

    print("Name1.2", "file_path", input1[2])
    file_path = input1[2]

    print("Name1.3", "sw_path", input1[3])
    sw_path = input1[3]

    print("Name2.0", "association_df", input2[0])
    association_df = pd.read_csv(
        input2[1], delimiter="\s+", dtype=str, header=None, names=["GO1", "GO2", "Ass"]
    )

    print("Name2.1", "category_df", input2[1])
    category_df = pd.read_csv(
        input2[2], delimiter="\s+", dtype=str, header=None, names=["GO", "Category"]
    )

    print("Name2.2", "child_parent_percentages_df", input2[2])
    child_parent_percentages_df = pd.read_csv(
        input2[3],
        delimiter="\s+",
        dtype=str,
        header=None,
        usecols=[1, 2, 3, 4, 5, 6, 7, 8],
        names=[
            "Child",
            "Parent",
            "Percentage",
            "None1",
            "None2",
            "None3",
            "None4",
            "None5",
        ],
    )

    print("Name2.3", "go_term1", input2[3])
    go_term1 = pd.read_csv(
        input2[4], delimiter="\t", dtype=str, header=None, names=["Id", "Name"]
    )

    print("Name2.4", "go_term2", input2[4])
    go_term2 = pd.read_csv(
        input2[4], delimiter="\t", dtype=str, header=None, names=["Id", "Name"]
    )

    print("Name2.5", "loc_2019_basic", input2[5])
    loc_2019_basic = input2[5]

    print("Name3", "df_len_seq", input3[0])
    df_len_seq = pd.read_csv(input3[0], dtype=str)

    print("Name3.1", "df_feat", input3[1])
    df_feat = pd.read_csv(input3[1], dtype=str, sep="\t")

    print("Name3.2", "base_url", input3[2])
    base_url = input3[2]

    output_dir_files = f"{org_var[0]}/{org_var[1]}"
    os.makedirs(output_dir_files, exist_ok=True)
    pygosemsim_handler = PyGOsemsimHandler()

    # Step1: Read in validation df
    val_or_test_df = pd.read_parquet(input1[0])
    val_or_test_df = (
        val_or_test_df[["GO_Term", "Ontology", "Raw_Score", "Rank", "Name"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    val_or_test_df = val_or_test_df.rename(columns={"GO_Term": "GO"}).reset_index(
        drop=True
    )
    single_predict_df, single_predict_prot_list = pygosemsim_handler.get_dataframe(
        val_or_test_df
    )
    del val_or_test_df
    print(
        "Step1 is a reminder that not all inputs are loaded in the begining of .py. Group 5 inputs are loaded right before they are used"
    )
    logger.info("Step1 is done")

    # Step2: Get the Ground Truth dataset, with 6 edited Species
    updated_uniprot_6_species_with_goterms_gaf = pd.read_csv(
        input1[1], sep=",", header=0, low_memory=False
    )
    updated_uniprot_6_species_with_g_oterms_gaf2 = (
        pygosemsim_handler.get_6_species_with_gaf(
            updated_uniprot_6_species_with_goterms_gaf
        )
    )
    print("Step2 is done")
    logger.info("Step2 is done")

    # Step3: Combine Entire SwissProt Dataset (this will be used to get all the GO terms and to propagate them to the root, in order to get
    # frequencies across the gene ontology. )
    simple_dic_var_partial1, final_go_partial_annot1, final_key_partial_annot1 = (
        pygosemsim_handler.combining_partial_files(file_path, sw_path)
    )
    print("Step3 is done")
    logger.info("Step3 is done")

    # Step4: I am simplifying the "simple_dic_var_partial1 dic". Simple_dic_var_partial2 is a same simple_dic_var_partial1 (keys are completly the same), except the structure is different (below)
    simple_dic_var_partial2 = defaultdict(lambda: defaultdict(dict))
    for (
        key,
        details,
    ) in (
        simple_dic_var_partial1.items()
    ):  # key is a protein and details are the go terms
        for go in details:
            simple_dic_var_partial2[key]["annotation"][go] = {go: None}
    print("Step4 is done")
    logger.info("Step4 is done")

    # Step5: I am filtering the 9006 proteins that I initially selected to correspond with the ~706 proteins that I run through the PFP and in val
    updated_uniprot_6_species_with_g_oterms_gaf2 = (
        updated_uniprot_6_species_with_g_oterms_gaf2[
            updated_uniprot_6_species_with_g_oterms_gaf2["Protein"].isin(
                single_predict_df["Protein"].drop_duplicates().to_list()
            )
        ]
    )
    print("Step5 is done")
    logger.info("Step5 is done")

    # Step6: I used pygosemsim's code, and I rewrote that to save the "Alt_Dic", "Global_Deleted" list, and dicts where the values were replaced
    G, global_alt_dic, global_repl_dic, global_ob_repl_dic, global_deleted = (
        graph.from_resource("OBO_2022_09_19_go-basic")
    )
    global_deleted1 = [el for el in global_deleted if el not in global_ob_repl_dic]
    print("Step6 is done")
    logger.info("Step6 is done")

    # Step7: Here, I am saving the 3 dictionaries, so I can always read them in
    with open(f"./{output_dir_files}/altD.pickle", "wb") as handle:
        pickle.dump(global_alt_dic, handle, protocol=4)
    with open(f"./{output_dir_files}/replD.pickle", "wb") as handle:
        pickle.dump(global_ob_repl_dic, handle, protocol=4)
    with open(f"./{output_dir_files}/deleted.pickle", "wb") as handle:
        pickle.dump(global_deleted1, handle, protocol=4)
    print("Step7 is done")
    logger.info("Step7 is done")

    # Step8: this function takes the GAF answers (for the training set) and propagates them with the "is_a" and "part_of" edges.
    # There is a similar fx to propagate the entire SwissProt
    updated_uniprot_6_species_with_goterms_gaf2_propagated = (
        pygosemsim_handler.propagate_until_root(
            G,
            global_alt_dic,
            global_repl_dic,
            global_deleted1,
            updated_uniprot_6_species_with_g_oterms_gaf2,
        )
    )
    print("Step8 is done")
    logger.info("Step8 is done")

    # Step9: Generated a Real and a Predicted GO terms answers
    real_go_dic1, pred_go_dic = (
        pygosemsim_handler.reshape_real_and_pred_propagated_go_terms(
            single_predict_df, updated_uniprot_6_species_with_goterms_gaf2_propagated
        )
    )
    print("Step9 is done")
    logger.info("Step9 is done")

    # Step10: I had to write 3 helper functions, because I wanted to exclude the GO terms that were not in the db from the correct answers, because of the slight time difference.
    list_combined_pfp = pygosemsim_handler.time_diff_check_1(
        association_df, category_df, child_parent_percentages_df, go_term1, go_term2
    )

    # Step10_b:find all the GO terms that are in the 2022 obo, but not in the files
    go_basic_obo_df_not_found_in2019 = pygosemsim_handler.time_diff_check_2(
        G, loc_2019_basic, list_combined_pfp, global_ob_repl_dic, global_alt_dic
    )

    # Step10_d:sort the real_go_dic1
    sorted_real_go_dic1 = pygosemsim_handler.sort_real_go_dict(real_go_dic1)

    # Step10_d:find all the GO terms that are in the 2022 obo, but not in the files and excluded those from the RealGO_dic1 (real answers)
    # I broke the fx into 2 parts: part1 and part2
    # Step10_d_1:call to the part1
    del_dict, add_dict = {}, {}
    sorted_real_go_dic1, return_del_dict, return_add_dict = (
        pygosemsim_handler.time_diff_check_3_update_real_go_dict1_part1(
            sorted_real_go_dic1, del_dict, add_dict, go_basic_obo_df_not_found_in2019
        )
    )
    # Step10_d_2:call to the part2
    real_go_dic1, list_removed_prop = (
        pygosemsim_handler.time_diff_check_3_update_real_go_dict1_part2(
            sorted_real_go_dic1, go_basic_obo_df_not_found_in2019
        )
    )  # actual call
    print("Step10 is done")
    logger.info("Step10 is done")

    # Step11: I propagated the SwissGO terms
    go_terms_propagation_dic_sw = pygosemsim_handler.propagate_gaf_until_root(
        simple_dic_var_partial2, G, global_alt_dic, global_ob_repl_dic, global_deleted1
    )
    print("Step11 is done")
    logger.info("Step11 is done")

    # Step12: This fx breaks the propagated SW into 3 ontologies
    joined_list_gaf_pb_sw, joined_list_gaf_cc_sw, joined_list_gaf_mf_sw = (
        pygosemsim_handler.breaking_propagation_sw_to_ont(go_terms_propagation_dic_sw)
    )
    print("Step12 is done")
    logger.info("Step12 is done")

    # Step13: This fx changes list of all GO terms per ontology to the df, and counts how many times each GO term is seen
    pb_sw_df = pygosemsim_handler.joined_lsit_fo_function(joined_list_gaf_pb_sw, "PB")
    cc_sw_df = pygosemsim_handler.joined_lsit_fo_function(joined_list_gaf_cc_sw, "CC")
    mf_sw_df = pygosemsim_handler.joined_lsit_fo_function(joined_list_gaf_mf_sw, "MF")
    print("Step13 is done")
    logger.info("Step13 is done")

    # Step14: Generating a dict, by comparing the df against the OBO (some GO terms may have count of 0 (they are not found in the df))
    dic_with_count_pb_sw = pygosemsim_handler.dic_with_count(G, pb_sw_df)
    dic_with_count_cc_sw = pygosemsim_handler.dic_with_count(G, cc_sw_df)
    dic_with_count_mf_sw = pygosemsim_handler.dic_with_count(G, mf_sw_df)
    print("Step14 is done")
    logger.info("Step14 is done")

    # Step15: This fx calculated the Information Concent (IC) by dividing how many times each GO term is seen against how many times root is seen
    # I am using Lin's score here (no need to normalize, also not choosing the biggest ancestor, but rather additive)
    dic2_with_i_bp_sw = pygosemsim_handler.get_ic(dic_with_count_pb_sw, "GO:0008150")
    dic2_with_i_cc_sw = pygosemsim_handler.get_ic(dic_with_count_cc_sw, "GO:0005575")
    dic2_with_i_mf_sw = pygosemsim_handler.get_ic(dic_with_count_mf_sw, "GO:0003674")

    # Step16_a: Just combine all 3 dictionaries
    dic_ic_bp_sw = dic2_with_i_bp_sw.copy()
    dic_ic_bp_sw.update(dic2_with_i_cc_sw)
    dic_ic_bp_sw.update(dic2_with_i_mf_sw)
    dic_ic_bp_cc_mf_sw = dic_ic_bp_sw.copy()
    print("Step15 is done")
    logger.info("Step15 is done")

    # Step16: Checkpoint: All of the sets should be empthy
    print(
        f"There are {pygosemsim_handler.checkpoint1(dic2_with_i_bp_sw, dic2_with_i_cc_sw, dic2_with_i_mf_sw)} intersecting GO terms within the 3 sets"
    )
    logger.info(
        f"There are {pygosemsim_handler.checkpoint1(dic2_with_i_bp_sw, dic2_with_i_cc_sw, dic2_with_i_mf_sw)} intersecting GO terms within the 3 sets"
    )
    print("Step16 is done")
    logger.info("Step16 is done")

    # Step 17: This is the big fx: For this fx, I am iterating over all the protein's, and generating a lin score for each prediction GO term against each real GO term, and I am choosing the biggest one of (Pred1 to Real1 is 1, pred1 to Real2 is 0.5: I would choose 1st relationship)
    # Step 17_a: Calling the main fx
    protein_ontologies = {}
    count = 0
    for prot in single_predict_prot_list[:]:
        count = count + 1
        (
            edited_temp_p_c,
            edited_temp_r_c,
            binary,
            real_missing,
            list_lin_tuples_c,
            list_lin_tuples_f,
            list_lin_tuples_p,
            temp_r,
            temp_r_c,
            temp_r_f,
            temp_r_p,
            temp_p,
            temp_p_c,
            temp_p_f,
            temp_p_p,
            prot,
        ) = pygosemsim_handler.get_ont_lists(
            prot,
            real_go_dic1,
            pred_go_dic,
            global_alt_dic,
            global_ob_repl_dic,
            global_deleted1,
        )
        protein_ontologies[prot] = {
            "binary": binary,
            "real_missing": real_missing,
            "list_lin_tuples_c": list_lin_tuples_c,
            "list_lin_tuples_f": list_lin_tuples_f,
            "list_lin_tuples_p": list_lin_tuples_p,
        }

    # Step 17_b: Just saving the results of the previous step
    with open(f"./{output_dir_files}/realGO_dic.pickle", "wb") as handle:
        pickle.dump(real_go_dic1, handle, protocol=4)
    with open(f"./{output_dir_files}/pred_go_dic.pickle", "wb") as handle:
        pickle.dump(pred_go_dic, handle, protocol=4)
    protein_ontologies1 = protein_ontologies.copy()
    with open(f"./{output_dir_files}/protein_ontologies1.pickle", "wb") as handle:
        pickle.dump(protein_ontologies1, handle, protocol=4)
    print("Step17 is done")
    logger.info("Step17 is done")

    # Step18: The "lin_score" is written to calculate the lin score and the "choosing_highest_score_within_pred_real" calculates the highval LinScore
    similarity.precalc_lower_bounds(G)
    dic_best_lin_score_per_go_c = (
        pygosemsim_handler.choosing_highest_score_within_pred_real(
            G,
            dic_ic_bp_cc_mf_sw,
            protein_ontologies1,
            "list_lin_tuples_c",
            go_basic_obo_df_not_found_in2019,
        )
    )
    dic_best_lin_score_per_go_f = (
        pygosemsim_handler.choosing_highest_score_within_pred_real(
            G,
            dic_ic_bp_cc_mf_sw,
            protein_ontologies1,
            "list_lin_tuples_f",
            go_basic_obo_df_not_found_in2019,
        )
    )
    dic_best_lin_score_per_goP = (
        pygosemsim_handler.choosing_highest_score_within_pred_real(
            G,
            dic_ic_bp_cc_mf_sw,
            protein_ontologies1,
            "list_lin_tuples_p",
            go_basic_obo_df_not_found_in2019,
        )
    )
    print("Step18 is done")
    logger.info("Step18 is done")

    # Step19: This fx changes a LinScore dict to the df
    lin_dataframe_c = pygosemsim_handler.lin_dic_to_df(dic_best_lin_score_per_go_c)
    lin_dataframe_f = pygosemsim_handler.lin_dic_to_df(dic_best_lin_score_per_go_f)
    lin_dataframe_p = pygosemsim_handler.lin_dic_to_df(dic_best_lin_score_per_goP)
    # display(lin_dataframe_p.tail(2))
    print("Ste19 is done")
    logger.info("Step19 is done")

    # Step20: This fx adds a root flag & combiens the 3 df
    lin_dataframe_c1 = pygosemsim_handler.add_root_flag(
        lin_dataframe_c, "GO:0005575", "C"
    )
    lin_dataframe_f1 = pygosemsim_handler.add_root_flag(
        lin_dataframe_f, "GO:0003674", "F"
    )
    lin_dataframe_p1 = pygosemsim_handler.add_root_flag(
        lin_dataframe_p, "GO:0008150", "P"
    )
    print("Step20 is done")
    logger.info("Step20 is done")

    # Step20_a:3 dfs are combined
    Lin_df_3 = pd.concat([lin_dataframe_c1, lin_dataframe_f1, lin_dataframe_p1], axis=0)
    Lin_df_3 = Lin_df_3.drop_duplicates()

    # Step21: Save the results
    with open(
        f"./{output_dir_files}/dic_best_lin_score_per_go_c.pickle", "wb"
    ) as handle:
        pickle.dump(dic_best_lin_score_per_go_c, handle, protocol=4)
    with open(
        f"./{output_dir_files}/dic_best_lin_score_per_go_f.pickle", "wb"
    ) as handle:
        pickle.dump(dic_best_lin_score_per_go_f, handle, protocol=4)
    with open(
        f"./{output_dir_files}/dic_best_lin_score_per_goP.pickle", "wb"
    ) as handle:
        pickle.dump(dic_best_lin_score_per_goP, handle, protocol=4)
    print("Step21 is done")
    logger.info("Step21 is done")

    # Step22: Pull out binary data, out of the protein_ontologies1 dict
    binary_dataframe = pygosemsim_handler.generate_binary_dataframe(protein_ontologies1)
    print("Step22 is done")
    logger.info("Step22 is done")

    # Step23_a: Checked the predicted df(from python4.pyPFP) against the mapped/replaced/deleted dictionaries
    single_predict_df1 = pygosemsim_handler.cleaning_predictions(
        single_predict_df, global_alt_dic, global_ob_repl_dic, global_deleted1
    )
    single_predict_df1.to_pickle(
        f"./{output_dir_files}/single_predict_df1.pkl", protocol=4
    )
    single_predict_df2 = single_predict_df1.drop(
        [
            "changed to replD",
            "changed to altD",
            "presence of deleted",
            "presence of replD",
            "presence of altD",
        ],
        axis=1,
    )

    # Step23_b:Merging the 2 df
    single_predict_df_b = pd.merge(
        single_predict_df2,
        binary_dataframe,
        left_on=["Protein", "changed to deleted"],
        right_on=["Protein", "PredictedGO"],
        how="inner",
    )
    single_predict_df_b = single_predict_df_b.drop(["PredictedGO"], axis=1)
    single_predict_df_b_l = pd.merge(
        single_predict_df_b,
        Lin_df_3,
        left_on=["Protein", "changed to deleted", "Predicted_Ontology"],
        right_on=["Protein", "PFP_Predicted_and_edited", "LinOnt1"],
        how="inner",
    )
    single_predict_df_b_l = single_predict_df_b_l.drop(
        ["PFP_Predicted_and_edited"], axis=1
    )
    single_predict_df_b_l = single_predict_df_b_l.drop(["LinOnt1"], axis=1)
    single_predict_df_b_l_sorted = single_predict_df_b_l.sort_values(
        by=["Protein", "PFP_Pred_GO_term", "Raw_Score", "Rank", "Raw_Score", "Rank"]
    ).reset_index(drop=True)
    print("Step23 is done")
    logger.info("Step23 is done")

    # Step24: Do merging/renaming operations and generate a list of missing feautures
    df_len_seq = df_len_seq[["Protein", "Seq_length"]]
    df_len_seq.rename(columns={"Protein": "Name"}, inplace=True)
    # annotation score
    df_feat.rename(columns={"Protein": "Name"}, inplace=True)
    protein_df = (
        single_predict_df_b_l_sorted[["Protein"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    protein_df_list = list(single_predict_df_b_l_sorted["Protein"].drop_duplicates())
    protein_df.rename(columns={"Protein": "Name"}, inplace=True)
    dfs = [protein_df, df_len_seq, df_feat]
    merged_df = reduce(
        lambda left, right: pd.merge(left, right, on="Name", how="inner"), dfs
    )
    merged_df_protein_list = list(merged_df["Name"].drop_duplicates())
    missing_len_feat = [x for x in protein_df_list if x not in merged_df_protein_list]
    print("Step24 is done")
    logger.info("Step24 is done")

    # Step25: Use that list of missing features, query the uniprot in Real Time, and again merge/rename/clean teh df
    result_df = pygosemsim_handler.query_uniprot_for_stats(base_url, missing_len_feat)
    feat_length_df = pd.concat([merged_df, result_df], axis=0)
    feat_length_df[["Seq_length", "AnnotationScore", "ProteinExistence"]] = (
        feat_length_df[["Seq_length", "AnnotationScore", "ProteinExistence"]].astype(
            float
        )
    )
    single_predict_df_b_l_sorted_f_l = pd.merge(
        single_predict_df_b_l_sorted,
        feat_length_df,
        left_on=["Protein"],
        right_on=["Name"],
        how="inner",
    )
    single_predict_df_b_l_sorted_f_l = single_predict_df_b_l_sorted_f_l
    single_predict_df_b_l_sorted_f_l["Predicted_Ontology_Retained"] = (
        single_predict_df_b_l_sorted_f_l["Predicted_Ontology"]
    )
    single_predict_df_b_l_sorted_f_l = pd.get_dummies(
        single_predict_df_b_l_sorted_f_l, columns=["Predicted_Ontology"]
    )
    single_predict_df_b_l_sorted_f_l = single_predict_df_b_l_sorted_f_l.copy(deep=True)
    print("Step25 is done")
    logger.info("Step25 is done")

    # Step26: Load the original PFP Prediciton's table as dict,
    print("Name5.1", "df_main_100", input5[2])
    # Comes from the Blast Data (Raw)
    df_main_100 = input5[2]
    training_handler = SparseTrainingArrayHandler()

    df_main_100 = training_handler.generate_original_predictions(df_main_100)
    df_main_100_df = pd.concat(df_main_100.values())
    df_main_100_df = df_main_100_df.drop("Verbal_Desc", axis=1)
    df_main_100_df["GO_Term"] = "GO:" + df_main_100_df["GO_Term"].astype(str)
    df_main_100_df["Ontology"] = df_main_100_df["Ontology"].str.upper()
    print("Step26 is done")
    logger.info("Step26 is done")

    # Step27: Process df: merge current df with the df_main_100, add depth and filter by dropping some columns
    single_predict_df_b_l_sorted_f_l1 = pd.merge(
        single_predict_df_b_l_sorted_f_l,
        df_main_100_df,
        left_on=[
            "Raw_Score",
            "Rank",
            "PFP_Pred_GO_term",
            "Predicted_Ontology_Retained",
            "Name",
        ],
        right_on=["Raw_Score", "Rank", "GO_Term", "Ontology", "Name"],
        how="inner",
    )
    single_predict_df_b_l_sorted_f_l1.rename(
        columns={
            "PFP_Pred_GO_term": "Original_PFP_Pred_GO_term",
            "changed to deleted": "Updated_PFP_Pred_GO_term",
        },
        inplace=True,
    )
    actual_list = single_predict_df_b_l_sorted_f_l1[
        "Updated_PFP_Pred_GO_term"
    ].to_list()
    actual_list2 = single_predict_df_b_l_sorted_f_l1["Predicted_Ontology_Retained"]
    single_predict_df_b_l_sorted_f_l_with_depth = pygosemsim_handler.add_go_term_depth(
        G, single_predict_df_b_l_sorted_f_l1
    )
    single_predict_df_b_l_sorted_f_l_with_depth = (
        single_predict_df_b_l_sorted_f_l_with_depth.drop(
            ["Predicted_Ontology_Retained", "Ontology", "GO_Term"], axis=1
        )
    )
    single_predict_df_b_l_sorted_f_l_with_depth.to_pickle(
        f"./{output_dir_files}/single_predict_df_b_l_sorted_f_l_with_depth.pkl",
        protocol=4,
    )
    print("Step27 is done")
    logger.info("Step27 is done")

    dataset_handler = TrainValidationSplitHandler()

    # Step28:Read in the initial dictionaries with the raw data
    all_dictionaries = dataset_handler.load_4_class_dictionaries(
        input6[0], input6[1], input6[2], input6[3]
    )
    # Step2:Separate main dict into 4 independent ones
    new_shape_before_hist_a = pd.concat(
        all_dictionaries["df_main_allGO_hist0"].values()
    )
    new_shape_before_hist_b = pd.concat(
        all_dictionaries["df_main_allAss_hist0"].values()
    )
    new_shape_before_hist_c = pd.concat(
        all_dictionaries["df_main_allCP_GO_hist0"].values()
    )
    new_shape_before_hist_d = pd.concat(
        all_dictionaries["df_main_allCP_Ass_hist0"].values()
    )
    print("Step28 is done")
    logger.info("Step28 is done")

    # Step29: Read the Validation Data in and Split the initial validation data by classes
    val_or_test_df = pd.read_parquet(input1[0])
    (
        pre_val_or_test_df_a,
        pre_val_or_test_df_b,
        pre_val_or_test_df_c,
        pre_val_or_test_df_d,
    ) = training_handler.split_dataframe_by_class(val_or_test_df)
    print("Step29 is done")
    logger.info("Step29 is done")

    # Step30: Using the original class-specific dataset, pick up class-specific features for the CV's dataset
    val_or_test_df_a = training_handler.enrich_dataframe_by_class(
        pre_val_or_test_df_a, new_shape_before_hist_a
    )
    val_or_test_df_b = training_handler.enrich_dataframe_by_class(
        pre_val_or_test_df_b, new_shape_before_hist_b
    )
    val_or_test_df_c = training_handler.enrich_dataframe_by_class(
        pre_val_or_test_df_c, new_shape_before_hist_c
    )
    val_or_test_df_d = training_handler.enrich_dataframe_by_class(
        pre_val_or_test_df_d, new_shape_before_hist_d
    )
    print("Step30 is done")
    logger.info("Step30 is done")

    # Step31: Load the "tmp4_to_analyze_df_train" and Add Raw Counts from the "val_or_test_df_a/B/C/D" to the train's generated bins
    # tmp4_to_analyze_df_train is a intermediate step in the training df (I am using it to get the bins's ranges)
    # Comes from the training data
    val_test_count_handler = ValTestCountHandler()
    tmp4_to_analyze_df_train = pd.read_pickle(input5[0])

    val_or_test_df_ind_classes_aa = val_test_count_handler.adding_raw_counts_to_bins(
        val_or_test_df_a, tmp4_to_analyze_df_train, "AA", "-log+b"
    )
    val_or_test_df_ind_classes_bb = val_test_count_handler.adding_raw_counts_to_bins(
        val_or_test_df_b, tmp4_to_analyze_df_train, "BB", "-log+b"
    )
    val_or_test_df_ind_classes_cc = val_test_count_handler.adding_raw_counts_to_bins(
        val_or_test_df_c, tmp4_to_analyze_df_train, "CC", "-log+b_P"
    )
    val_or_test_df_ind_classes_dd = val_test_count_handler.adding_raw_counts_to_bins(
        val_or_test_df_d, tmp4_to_analyze_df_train, "DD", "-log+b_P"
    )

    print("Step31 is done")
    logger.info("Step31 is done")

    # Step32: Combining Val df of 4 diff classes with each other by "Outer Merge" **Normally it is inner
    tmp1 = pd.merge(
        val_or_test_df_ind_classes_aa,
        val_or_test_df_ind_classes_bb,
        on=["GO_Term", "Ontology", "Raw_Score", "Rank", "Name"],
        how="outer",
    )
    tmp1 = pd.merge(
        tmp1,
        val_or_test_df_ind_classes_cc,
        on=["GO_Term", "Ontology", "Raw_Score", "Rank", "Name"],
        how="outer",
    )
    tmp1 = pd.merge(
        tmp1,
        val_or_test_df_ind_classes_dd,
        on=["GO_Term", "Ontology", "Raw_Score", "Rank", "Name"],
        how="outer",
    )
    tmp1["GO_Term"] = "GO:" + tmp1["GO_Term"]
    tmp1["Ontology"] = tmp1["Ontology"].str.upper()
    tmp1 = tmp1.fillna(0)
    print("Step32 is done")
    logger.info("Step32 is done")

    # Step33: Merging initial df with the tmp1 (created in the previous step with the counts data), and cleaning some rows
    single_predict_df_b_l_sorted_f_l_with_depth2 = pd.merge(
        tmp1,
        single_predict_df_b_l_sorted_f_l_with_depth,
        left_on=["GO_Term", "Raw_Score", "Rank", "Name"],
        right_on=["Original_PFP_Pred_GO_term", "Raw_Score", "Rank", "Name"],
        how="inner",
    )
    if list(single_predict_df_b_l_sorted_f_l_with_depth2["Name"]) == list(
        single_predict_df_b_l_sorted_f_l_with_depth2["Protein"]
    ):
        single_predict_df_b_l_sorted_f_l_with_depth2 = (
            single_predict_df_b_l_sorted_f_l_with_depth2.drop("Protein", axis=1)
        )
    if list(single_predict_df_b_l_sorted_f_l_with_depth2["GO_Term"]) == list(
        single_predict_df_b_l_sorted_f_l_with_depth2["Original_PFP_Pred_GO_term"]
    ):
        single_predict_df_b_l_sorted_f_l_with_depth2 = (
            single_predict_df_b_l_sorted_f_l_with_depth2.drop("GO_Term", axis=1)
        )
    pre_final_val_or_test_df = single_predict_df_b_l_sorted_f_l_with_depth2.sort_values(
        "Raw_Score"
    ).reset_index(drop=True)

    pre_final_val_or_test_df.to_pickle(
        f"./{output_dir_files}/pre_final_val_or_test_df.pkl", protocol=4
    )
    print("Step33 is done")
    logger.info("Step33 is done")

    # Step34: Load the "final_train_df" and Get Rid of the Roots, changing the column's name and reorganizing the val df to look like training df
    # The reference df used in the fx is the final_training_df, but I am only using it for the col's names and col's order. No leakage!
    # Comes from the training data

    final_train_df = pd.read_pickle(input5[1])
    final_val_or_test_df = val_test_count_handler.clean_val_or_df(
        pre_final_val_or_test_df,
        final_train_df,
        ["Original_PFP_Pred_GO_term", "Updated_PFP_Pred_GO_term"],
        ["GO:0008150", "GO:0005575", "GO:0003674"],
    )

    final_val_or_test_df.to_pickle(
        f"./{output_dir_files}/final_{org_var[4]}_df.pkl", protocol=4
    )
    print("Step34 is done")
    logger.info("Step34 is done")
    print("Finished the main script", flush=True)
    logger.info("Global1_part456.py is finished")
    gc.collect()
    logging.shutdown()


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(
            "Usage: python script.py <input1> <input2> <input3>, <org_var>, <input5>, <input6>"
        )
        sys.exit(1)

    input1_arg = sys.argv[
        1
    ]  # Collect all arguments after the script name. Those arguments are the original 4 df
    input1 = [col.strip() for col in input1_arg.split(",")]
    input2_arg = sys.argv[
        2
    ]  # Collect all arguments after the script name. Those arguments are needed
    input2 = [col.strip() for col in input2_arg.split(",")]
    input3_arg = sys.argv[
        3
    ]  # Collect all arguments after the script name. Those arguments are needed
    input3 = [col.strip() for col in input3_arg.split(",")]
    input4_arg = sys.argv[
        4
    ]  # Collect all arguments after the script name. Those arguments are the CV pair
    org_var = [col.strip() for col in input4_arg.split(",")]
    input5_arg = sys.argv[
        5
    ]  # Collect all arguments after the script name. Those arguments are the needed to send the raw counts into the training bins, and reformat final_val_or_test_df's columns to be in the same order as final_train_df
    input5 = [col.strip() for col in input5_arg.split(",")]
    input6_arg = sys.argv[
        6
    ]  # Collect all arguments after the script name. Those arguments are the needed to send the raw counts into the training bins, and reformat final_val_or_test_df's columns to be in the same order as final_train_df
    input6 = [col.strip() for col in input6_arg.split(",")]
    mainprocess(input1, input2, input3, org_var, input5, input6)
