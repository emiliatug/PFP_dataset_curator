import gc
import pickle
from functools import reduce

import pandas as pd

# pd.set_option('future.no_silent_downcasting', True)
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)

from Modules.pygosemsim_handler import *
from Utils.utils import *
from pygosemsim import graph
from pygosemsim import similarity


def mainprocess(
    input1: List[str], input2: List[str], input3: List[str], org_var: List[str]
) -> None:
    """
    Main processing function for Global Part1: Pygosemsim Integrator.

    Parameters:
    single_prediction_inputs (List[str]): List of paths related to single prediction inputs.
    association_category_inputs (List[str]): List of paths for association and category data.
    sequence_feature_inputs (List[str]): List of paths for sequence length and feature data.
    run_metadata (List[str]): List of metadata including output directories and run identifiers.
    """
    print("output_dir", org_var)
    output_dir = org_var[2]
    os.makedirs(output_dir, exist_ok=True)
    print("output_dir", output_dir)

    # Part0.5 Set up logger file
    logger = generate_logging_file(
        name=f"{org_var[3]}_{org_var[1]}", loc_to_write=output_dir
    )
    # logger=generate_logging_file(name=f'{org_var[3]_{org_var[1]}}', loc_to_write=output_dir)
    logger.info("Starting the Global Part1: Data Generation:Train: Part5")

    # Step0: Real in all the inputs
    print("All inputs are: ")
    print("Name1.0", "single_predict_df", input1[0])
    single_predict_df = pd.read_pickle(input1[0])

    print("Name1.1", "updated_uniprot_6_species_with_GOterms_GAF", input1[1])

    print("Name1.2", "file_path", input1[2])
    file_path = input1[2]

    print("Name1.3", "sw_path", input1[3])
    sw_path = input1[3]

    print("Name2.0", "association_df", input2[0])
    association_df = pd.read_csv(
        input2[0], delimiter="\s+", dtype=str, header=None, names=["GO1", "GO2", "Ass"]
    )

    print("Name2.1", "category_df", input2[1])
    category_df = pd.read_csv(
        input2[1], delimiter="\s+", dtype=str, header=None, names=["GO", "Category"]
    )

    print("Name2.2", "child_parent_percentages_df", input2[2])
    child_parent_percentages_df = pd.read_csv(
        input2[2],
        delimiter="\s+",
        dtype=str,
        header=None,
        usecols=[0, 1, 2, 3, 4, 5, 6, 7],
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

    print("Name2.3", "go_term_1", input2[3])
    go_term_1 = pd.read_csv(
        input2[3], delimiter="\t", dtype=str, header=None, names=["Id", "Name"]
    )

    print("Name2.4", "go_term_2", input2[4])
    go_term_2 = pd.read_csv(
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

    handler = PyGOsemsimHandler()

    # Step1: Read in results of 4th step
    from_part4single_predict_df = single_predict_df.reset_index(drop=True)
    from_part4single_predict_df, from_part4single_predict_prot_list = (
        handler.get_dataframe(from_part4single_predict_df)
    )
    print("Step1 is done")
    logger.info("Step1 is done")

    # Step2: Get the Ground Truth dataset, with 6 edited Species
    updated_uniprot_6_species_with_GOterms_GAF = pd.read_csv(
        input1[1], sep=",", header=0, low_memory=False
    )
    updated_uniprot_6_species_with_GOterms_GAF2 = handler.get_6_species_with_gaf(
        updated_uniprot_6_species_with_GOterms_GAF
    )
    print("Step2 is done")
    logger.info("Step2 is done")

    # Step3: Combine Entire SwissProt Dataset (this will be used to get all the GO terms and to propagate them to the root, in order to get
    # frequencies across the gene ontology. )
    simple_dic_var_partial1, final_go_partial_annot1, final_key_partial_annot1 = (
        handler.combining_partial_files(file_path, sw_path)
    )
    print("Step3 is done")
    logger.info("Step3 is done")

    # Step4: I am simplying the "simple_dic_var_partial1 dic". Simple_dic_var_partial2 is a same simple_dic_var_partial1 (keys are completly the same), except the structure is different (below)
    simple_dic_var_partial2 = defaultdict(lambda: defaultdict(dict))
    for (
        key,
        details,
    ) in (
        simple_dic_var_partial1.items()
    ):  # key is a protein andn details are the go terms
        for go in details:
            simple_dic_var_partial2[key]["annotation"][go] = {go: None}
    print("Step4 is done")
    logger.info("Step4 is done")

    # Step5: Filtering the 9006 proteins that I initially selected to correspond with the ~3800 proteins that I run through the PFP
    updated_uniprot_6_species_with_GOterms_GAF2 = (
        updated_uniprot_6_species_with_GOterms_GAF2[
            updated_uniprot_6_species_with_GOterms_GAF2["Protein"].isin(
                from_part4single_predict_df["Protein"].drop_duplicates().to_list()
            )
        ]
    )
    print("Step5 is done")
    logger.info("Step5 is done")

    # Step6: Used pygosemsim's code, and restructured that to save the "Alt_Dic", "Global_Deleted" list, and dicts where the values were replaced
    # This should print ("format-version: 1.2","final version","43335","['part_of', 'is_a']")
    # G is a graph
    G, global_alt_dic, global_repl_dic, global_ob_repl_dic, global_deleted = (
        graph.from_resource("OBO_2022_09_19_go-basic")
    )
    global_deleted1 = [el for el in global_deleted if el not in global_ob_repl_dic]
    print("Step6 is done")
    logger.info("Step6 is done")

    # Step7: Saving the 3 dictionaries.
    with open(f"./{output_dir_files}/altD.pickle", "wb") as f:
        pickle.dump(global_alt_dic, f, protocol=4)
    with open(f"./{output_dir_files}/replD.pickle", "wb") as f:
        pickle.dump(global_ob_repl_dic, f, protocol=4)
    with open(f"./{output_dir_files}/deleted.pickle", "wb") as f:
        pickle.dump(global_deleted1, f, protocol=4)
    print("Step7 is done")
    logger.info("Step7 is done")

    # Step8: this function takes the GAF answers (for the training set) and propagates them with the "is_a" and "part_of" edges.
    # There is a similar fx to propagate the entire SwissProt
    updated_uniprot_6_species_with_GOterms_GAF2_propogated = (
        handler.propogate_until_root(
            G,
            global_alt_dic,
            global_repl_dic,
            global_deleted1,
            updated_uniprot_6_species_with_GOterms_GAF2,
        )
    )
    print("Step8 is done")
    logger.info("Step8 is done")

    # Step9: Generated a Real and a Predicted GO terms answers
    real_go_dic1, pred_go_dic = handler.reshape_real_and_pred_propagated_go_terms(
        from_part4single_predict_df,
        updated_uniprot_6_species_with_GOterms_GAF2_propogated,
    )
    print("Step9 is done")
    logger.info("Step9 is done")

    # Step10: I had to write 3 helper functions, because I wanted to exclude the GO terms that were not in the db from the correct answers, because of the slight time difference.
    # Step10_a:find all the GO terms that were occurring in the files
    list_combined_pfp = handler.time_diff_check_1(
        association_df, category_df, child_parent_percentages_df, go_term_1, go_term_2
    )

    # Step10_b:find all the GO terms that are in the 2022 obo, but not in the files
    go_basic_obo_df_not_found_in2019 = handler.time_diff_check_2(
        G, loc_2019_basic, list_combined_pfp, global_ob_repl_dic, global_alt_dic
    )

    # Step10_d:sort the real_go_dic1
    sorted_realGO_dic1 = handler.sort_real_go_dict(real_go_dic1)

    # Step10_d:find all the GO terms that are in the 2022 obo, but not in the files and excluded those from the RealGO_dic1 (real answers)
    # I broke the fx into 2 parts: part1 and part2
    # Step10_d_1:call to the part1
    del_dict, add_dict = {}, {}
    sorted_realGO_dic1, return_del_dict, return_add_dict = (
        handler.time_diff_check_3_update_real_go_dict1_part1(
            sorted_realGO_dic1, del_dict, add_dict, go_basic_obo_df_not_found_in2019
        )
    )
    # Step10_d_2:call to the part2
    real_go_dic1, list_removed_prop = (
        handler.time_diff_check_3_update_real_go_dict1_part2(
            sorted_realGO_dic1, go_basic_obo_df_not_found_in2019
        )
    )  # actual call
    print("Step10 is done")
    logger.info("Step10 is done")

    # Step11: I propagated the SwissGO terms
    go_terms_propagation_dic_sw = handler.propogate_gaf_until_root(
        simple_dic_var_partial2, G, global_alt_dic, global_ob_repl_dic, global_deleted1
    )
    print("Step11 is done")
    logger.info("Step11 is done")

    # Step12: This fx breaks the propagated SW into 3 ontologies
    joined_list_gaf_pb_sw, joined_list_gaf_cc_sw, joined_list_gaf_mf_sw = (
        handler.breaking_propog_sw_to_ont(go_terms_propagation_dic_sw)
    )
    print("Step12 is done")
    logger.info("Step12 is done")

    # Step13: This fx changes list of all GO terms per ontology to the df, and counts how many times each GO term is seen
    pb_sw_df = handler.joined_lsit_fo_function(joined_list_gaf_pb_sw, "PB")
    cc_sw_df = handler.joined_lsit_fo_function(joined_list_gaf_cc_sw, "CC")
    mf_sw_df = handler.joined_lsit_fo_function(joined_list_gaf_mf_sw, "MF")
    print("Step13 is done")
    logger.info("Step13 is done")

    # Step14: Generating a dict, by comparing the df against the OBO (some GO terms may have count of 0 (they are not found in the df))
    dic_with_count_pb_sw = handler.dic_with_count(G, pb_sw_df)
    dic_with_count_cc_sw = handler.dic_with_count(G, cc_sw_df)
    dic_with_count_mf_sw = handler.dic_with_count(G, mf_sw_df)
    print("Step14 is done")
    logger.info("Step14 is done")

    # Step15: This fx calculated the Information Concent (IC) by dividing how many times each GO term is seen against how many times root is seen
    # I am using Lin's score here (no need to normalize, also not choosing the biggest ancestor, but rather additive)
    dic2_with_i_bp_sw = handler.get_ic(dic_with_count_pb_sw, "GO:0008150")
    dic2_with_i_cc_sw = handler.get_ic(dic_with_count_cc_sw, "GO:0005575")
    dic2_with_i_mf_sw = handler.get_ic(dic_with_count_mf_sw, "GO:0003674")

    # Step16_a: Just combine all 3 dictionaries
    dic_ic_bp_sw = dic2_with_i_bp_sw.copy()
    dic_ic_bp_sw.update(dic2_with_i_cc_sw)
    dic_ic_bp_sw.update(dic2_with_i_mf_sw)
    dic_ic_bp_cc_mf_sw = dic_ic_bp_sw.copy()
    print("Step15 is done")
    logger.info("Step15 is done")

    # Step16: Checkpoint: All of the sets should be empty
    print(
        f"There are {handler.checkpoint1(dic2_with_i_bp_sw, dic2_with_i_cc_sw, dic2_with_i_mf_sw)} intersecting GO terms within the 3 sets"
    )
    logger.info(
        f"There are {handler.checkpoint1(dic2_with_i_bp_sw, dic2_with_i_cc_sw, dic2_with_i_mf_sw)} intersecting GO terms within the 3 sets"
    )
    print("Step16 is done")
    logger.info("Step16 is done")

    # Step 17: This is the big fx: For this fx, I am iterating over all the protein's, and generating a lin score for each prediction GO term against each real GO term, and I am choosing the biggest one of (Pred1 to Real1 is 1, pred1 to Real2 is 0.5: I would choose 1st relationship)
    # Step 17_a: Calling the main fx
    protein_ontologies = {}
    count = 0
    for prot in from_part4single_predict_prot_list[:]:
        # print(prot)
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
        ) = handler.get_ont_lists(
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
    with open(f"./{output_dir_files}/realGO_dic.pickle", "wb") as f:
        pickle.dump(real_go_dic1, f, protocol=4)
    with open(f"./{output_dir_files}/pred_go_dic.pickle", "wb") as f:
        pickle.dump(pred_go_dic, f, protocol=4)
    protein_ontologies1 = protein_ontologies.copy()
    with open(f"./{output_dir_files}/protein_ontologies1.pickle", "wb") as f:
        pickle.dump(protein_ontologies1, f, protocol=4)
    print("Step17 is done")
    logger.info("Step17 is done")

    # Step19: The "lin_score" is written to calculate the lin score and the "choosing_highest_score_within_pred_real" calculates the hightest LinScore
    # from pygosemsim import similarity
    similarity.precalc_lower_bounds(G)
    print("pass1 after17")
    dic_best_lin_score_per_go_c = handler.choosing_highest_score_within_pred_real(
        G,
        dic_ic_bp_cc_mf_sw,
        protein_ontologies1,
        "list_lin_tuples_c",
        go_basic_obo_df_not_found_in2019,
    )
    dic_best_lin_score_per_go_f = handler.choosing_highest_score_within_pred_real(
        G,
        dic_ic_bp_cc_mf_sw,
        protein_ontologies1,
        "list_lin_tuples_f",
        go_basic_obo_df_not_found_in2019,
    )
    dic_best_lin_score_per_go_p = handler.choosing_highest_score_within_pred_real(
        G,
        dic_ic_bp_cc_mf_sw,
        protein_ontologies1,
        "list_lin_tuples_p",
        go_basic_obo_df_not_found_in2019,
    )
    print("Step18 is done")
    logger.info("Step18 is done")

    # Step19: This fx changes a LinScore dict to the df
    lin_dataframe_c = handler.lin_dic_to_df(dic_best_lin_score_per_go_c)
    lin_dataframe_f = handler.lin_dic_to_df(dic_best_lin_score_per_go_f)
    lin_dataframe_p = handler.lin_dic_to_df(dic_best_lin_score_per_go_p)
    # display(lin_dataframe_p.tail(2))
    print("Ste19 is done")
    logger.info("Step19 is done")

    # Step20: This fx adds a root flag & combined the 3 df
    lin_dataframe_c1 = handler.add_root_flag(lin_dataframe_c, "GO:0005575", "C")
    lin_dataframe_f1 = handler.add_root_flag(lin_dataframe_f, "GO:0003674", "F")
    lin_dataframe_p1 = handler.add_root_flag(lin_dataframe_p, "GO:0008150", "P")

    # Step20_a:3 dfs are combined
    lin_df_3 = pd.concat([lin_dataframe_c1, lin_dataframe_f1, lin_dataframe_p1], axis=0)
    lin_df_3 = lin_df_3.drop_duplicates()
    # display(lin_df_3.head(2))
    print("Step20 is done")
    logger.info("Step20 is done")

    # Step21: Save the results
    with open(f"./{output_dir_files}/dic_best_lin_score_per_go_c.pickle", "wb") as f:
        pickle.dump(dic_best_lin_score_per_go_c, f, protocol=4)
    with open(f"./{output_dir_files}/dic_best_lin_score_per_go_f.pickle", "wb") as f:
        pickle.dump(dic_best_lin_score_per_go_f, f, protocol=4)
    with open(f"./{output_dir_files}/dic_best_lin_score_per_go_p.pickle", "wb") as f:
        pickle.dump(dic_best_lin_score_per_go_p, f, protocol=4)
    print("Step21 is done")
    logger.info("Step21 is done")

    # Step22: Pull out binary data, out of the protein_ontologies1 dict
    binary_dataframe = handler.generate_binary_dataframe(protein_ontologies1)
    print("Step22 is done")
    logger.info("Step22 is done")

    # Step23_a: Checked the predicted df(from python4.pyPFP) against the mapped/raplaced/deleted dictionaries
    from_part4single_predict_df1 = handler.cleaning_predictions(
        from_part4single_predict_df, global_alt_dic, global_ob_repl_dic, global_deleted1
    )

    from_part4single_predict_df1.to_pickle(
        f"./{output_dir_files}/from_part4single_predict_df1.pkl", protocol=4
    )
    from_part4single_predic_df2 = from_part4single_predict_df1.drop(
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
    from_part4single_predict_df_b = pd.merge(
        from_part4single_predic_df2,
        binary_dataframe,
        left_on=["Protein", "changed to deleted"],
        right_on=["Protein", "PredictedGO"],
        how="inner",
    )
    from_part4single_predict_df_b = from_part4single_predict_df_b.drop(
        ["PredictedGO"], axis=1
    )
    from_part4single_predict_df_b_merged = pd.merge(
        from_part4single_predict_df_b,
        lin_df_3,
        left_on=["Protein", "changed to deleted", "Predicted_Ontology"],
        right_on=["Protein", "PFP_Predicted_and_edited", "LinOnt1"],
        how="inner",
    )
    from_part4single_predict_df_b_merged = from_part4single_predict_df_b_merged.drop(
        ["PFP_Predicted_and_edited"], axis=1
    )
    from_part4single_predict_df_b_merged = from_part4single_predict_df_b_merged.drop(
        ["LinOnt1"], axis=1
    )
    from_part4single_predict_df_b_merged_sorted = (
        from_part4single_predict_df_b_merged.sort_values(
            by=["Protein", "PFP_Pred_GO_term", "Raw_Score", "Rank", "Raw_Score", "Rank"]
        ).reset_index(drop=True)
    )
    print("Step23 is done")
    logger.info("Step23 is done")

    # Step24: Do merging/renaming operations and generate a list of missing features
    df_len_seq = df_len_seq[["Protein", "Seq_length"]]
    df_len_seq.rename(columns={"Protein": "Name"}, inplace=True)
    # annotation score
    df_feat.rename(columns={"Protein": "Name"}, inplace=True)

    protein_df = (
        from_part4single_predict_df_b_merged_sorted[["Protein"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    protein_df_list = list(
        from_part4single_predict_df_b_merged_sorted["Protein"].drop_duplicates()
    )
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
    result_df = handler.query_uniprot_for_stats(base_url, missing_len_feat)
    feat_length_df = pd.concat([merged_df, result_df], axis=0)
    feat_length_df[["Seq_length", "AnnotationScore", "ProteinExistence"]] = (
        feat_length_df[["Seq_length", "AnnotationScore", "ProteinExistence"]].astype(
            float
        )
    )
    from_part4single_predict_df_b_l_sorted_f_l = pd.merge(
        from_part4single_predict_df_b_merged_sorted,
        feat_length_df,
        left_on=["Protein"],
        right_on=["Name"],
        how="inner",
    )
    from_part4single_predict_df_b_l_sorted_f_l = (
        from_part4single_predict_df_b_l_sorted_f_l
    )
    from_part4single_predict_df_b_l_sorted_f_l["Predicted_Ontology_Retained"] = (
        from_part4single_predict_df_b_l_sorted_f_l["Predicted_Ontology"]
    )
    from_part4single_predict_df_b_l_sorted_f_l = pd.get_dummies(
        from_part4single_predict_df_b_l_sorted_f_l, columns=["Predicted_Ontology"]
    )
    from_part4single_predict_df_b_l_sorted_f_l = (
        from_part4single_predict_df_b_l_sorted_f_l.copy(deep=True)
    )
    logger.info("Step25 is done")

    # Step26: Merge/sort/rename/and save
    from_part4single_predict_df_b_l_sorted_f_l.rename(
        columns={
            "PFP_Pred_GO_term": "Original_PFP_Pred_GO_term",
            "changed to deleted": "Updated_PFP_Pred_GO_term",
        },
        inplace=True,
    )
    from_part4single_predict_df_b_l_sorted_f_l_with_depth = handler.add_go_term_depth(
        G, from_part4single_predict_df_b_l_sorted_f_l
    )
    from_part4single_predict_df_b_l_sorted_f_l_with_depth.to_pickle(
        f"./{output_dir_files}/from_part4single_predict_df_b_l_sorted_f_l_with_depth.pkl",
        protocol=4,
    )
    print("Step26 is done")
    print("Finished the main script", flush=True)
    logger.info("Step26 is done")
    logger.info("Global1_part5.py is finished")
    gc.collect()
    logging.shutdown()


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input1> <input2> <input3>,<org_var>")
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
    mainprocess(input1, input2, input3, org_var)
