# import gc

from Modules.initial_preprocessing import InitialPreprocessor
from Utils.utils import *


def mainprocess(
    original_databases_dir: List[str],
    uniprot_url: str,
    intermediate_dict_dir: List[str],
    child_parent_dir: List[str],
    four_classes_dir: List[str],
    output_log_vars: List[str],
) -> None:
    output_dir = output_log_vars[0]
    os.makedirs(output_dir, exist_ok=True)

    logger = generate_logging_file(
        name=f"{output_log_vars[1]}", loc_to_write=output_dir
    )
    logger.info("Starting the Class Generation: Part1")

    preprocessor = InitialPreprocessor()

    # Step1: Download PFP Data
    df_main_100 = preprocessor.load_pfp_predictions(original_databases_dir[0])
    logger.info("Step 1 is done")

    # Step2: Download Data with the E values
    df_prot_eval_go = preprocessor.load_protein_evalues(
        original_databases_dir[1], df_main_100
    )
    logger.info("Step 2 is done")

    # Step3: Get the verbal name for the proteins
    df_names, protein_list = preprocessor.extract_verbal_names(
        original_databases_dir[2], df_prot_eval_go
    )
    logger.info("Step 3 is done")

    # Step4: Query the uniprot to get the missing verbal names, for some of the proteins
    df_names1 = preprocessor.query_uniprot_for_names(
        protein_list, uniprot_url, df_names
    )

    # Step5: Adding the verbal names derived from the 2 diff sources, and splitting a single dict to 3 diff dict
    df_prot_eval_go_upd, df_prot_eval_go_upd_go, df_prot_eval_go_upd_ass = (
        preprocessor.enrich_with_verbal_names(df_names1, df_prot_eval_go)
    )
    logger.info("Step 5 is done")

    saving_lists_of_locats_and_dicts(
        [intermediate_dict_dir[0], intermediate_dict_dir[1], intermediate_dict_dir[2]],
        [df_prot_eval_go_upd, df_prot_eval_go_upd_go, df_prot_eval_go_upd_ass],
    )
    logger.info("Step 6 is done")

    # Step7: This fx generates first 2 of 4 classes of dict
    df_main_all_go = preprocessor.merge_with_pfp(df_prot_eval_go_upd_go, df_main_100)
    df_main_all_ass = preprocessor.merge_with_pfp(df_prot_eval_go_upd_ass, df_main_100)
    logger.info("Step 7 is done")

    # Step8: generates dict of child-parent GO terms
    dict_ch_par = preprocessor.load_child_parent_relations(
        child_parent_dir[0], df_main_100
    )
    logger.info("Step 8 is done")

    # Step9: Generates class C and D
    df_main_all_cp_go = preprocessor.enrich_with_parental_relations(
        dict_ch_par, df_prot_eval_go_upd_go, df_main_100
    )
    df_main_all_cp_ass = preprocessor.enrich_with_parental_relations(
        dict_ch_par, df_prot_eval_go_upd_ass, df_main_100
    )
    logger.info("Step 9 is done")

    saving_lists_of_locats_and_dicts(
        [
            four_classes_dir[0],
            four_classes_dir[1],
            four_classes_dir[2],
            four_classes_dir[3],
        ],
        [df_main_all_go, df_main_all_ass, df_main_all_cp_go, df_main_all_cp_ass],
    )
    logger.info("Step 10 is done")

    # Step11: reorganize child-parent percentages
    preprocessor.organize_child_parent_percentages(
        child_parent_dir[1], child_parent_dir[2]
    )
    logger.info("Step 11 is done")

    # Step12: Check BLASTvsAss_Flag
    preprocessor.validate_single_flag(df_main_all_go, "BLASTvsAss_Flag", "True")
    preprocessor.validate_single_flag(df_main_all_go, "BLASTvsAss_Flag", "False")
    preprocessor.validate_single_flag(df_main_all_go, "BLASTvsAss_Flag", "True")
    preprocessor.validate_single_flag(df_main_all_go, "BLASTvsAss_Flag", "False")
    logger.info("Checkpoint1: all classes have appropriate BLASTvsAss_Flag")
    logger.info("Step 12 is done")

    # Step13: Check if all proteins from PFP exist in merged sets
    preprocessor.checkpoint_results(
        df_main_100,
        df_main_all_cp_go,
        df_main_all_cp_ass,
        df_main_all_go,
        df_main_all_ass,
    )
    logger.info(
        "Checkpoint2: All original results from the PFP are retrieved and match exactly within the 4 classes"
    )
    logger.info("Step 13 is done")

    logger.info("class_generator.py is done")


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(
            "Usage: python script.py <original_databases>, <uniprot_url>, <intermediate_dict_dir>, <child_parent_dir>,<four_classes_dir>,<output_log_vars>"
        )
        sys.exit(1)

    original_databases = [col.strip() for col in sys.argv[1].split(",")]
    uniprot_url = sys.argv[2]
    intermediate_dict_dir = [col.strip() for col in sys.argv[3].split(",")]
    child_parent_dir = [col.strip() for col in sys.argv[4].split(",")]
    four_classes_dir = [col.strip() for col in sys.argv[5].split(",")]
    output_log_vars = [col.strip() for col in sys.argv[6].split(",")]
    mainprocess(
        original_databases,
        uniprot_url,
        intermediate_dict_dir,
        child_parent_dir,
        four_classes_dir,
        output_log_vars,
    )
