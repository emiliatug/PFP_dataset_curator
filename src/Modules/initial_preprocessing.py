import gc
import glob
from typing import Dict, Any, Set, Tuple

import numpy as np
import requests

from Utils.utils import *


class InitialPreprocessor:
    """
    Object-oriented implementation of initial_preprocessing.py
    """

    def __init__(self):
        pass

    @staticmethod
    def load_pfp_predictions(input_dir: str) -> Dict[str, pd.DataFrame]:
        """
        Loads original PFP predictions.
        """
        df_main_100 = {}
        for filename in os.listdir(input_dir):
            if filename.endswith(".csv"):
                file_path = os.path.join(input_dir, filename)
                df = pd.read_csv(
                    file_path,
                    sep="\t",
                    dtype=str,
                    names=[
                        "GO_Term",
                        "Ontology",
                        "Out_of_Total",
                        "Raw_Score",
                        "Rank",
                        "Verbal_Desc",
                    ],
                )
                df["Raw_Score"] = pd.to_numeric(df["Raw_Score"])
                df_filtered = df[df["Raw_Score"] > 100.00]
                if not df_filtered.empty:
                    protein_name = filename[:-4]
                    df_filtered["Name"] = protein_name
                    df_main_100[protein_name] = df_filtered
        return df_main_100

    @staticmethod
    def load_protein_evalues(
        input_dir: str, df_main_100: Dict[str, pd.DataFrame]
    ) -> Dict[str, pd.DataFrame]:
        """
        Reads ProteinEvaluesGOterms data.
        """
        df_ProtEvalGO = {}
        for filename in os.listdir(input_dir):
            name = filename[:-11]
            if name in df_main_100.keys() and filename.endswith(".csv"):
                file_path = os.path.join(input_dir, filename)
                df = pd.read_csv(
                    file_path,
                    dtype=str,
                    sep=",",
                    names=[
                        "Type_of_GO",
                        "Actual_Blasted_GO_Term",
                        "found_protein",
                        "found_evalue",
                        "found_or_ass_go",
                    ],
                )
                df["BLASTvsAss_Flag"] = df["Type_of_GO"] == "TypeisGO"
                df["found_or_ass_go"] = df["found_or_ass_go"].str[:-1]
                df["Name"] = name
                df_ProtEvalGO[name] = df
        return df_ProtEvalGO

    @staticmethod
    def extract_verbal_names(
        base_pattern: str, df_prot_eval_go: Dict[str, pd.DataFrame]
    ) -> Tuple[pd.DataFrame, Set[str]]:
        """
        Extracts verbal names from BLAST files.
        """
        base_pattern = os.path.expanduser(base_pattern)
        data_dict = {}
        data_list = []
        for base_dir in glob.glob(base_pattern):
            file_pattern = os.path.join(base_dir, "*", "blast", "*_blast")
            for file_path in glob.glob(file_pattern):
                with open(file_path, "r") as file:
                    for line in file:
                        if line.startswith("sp"):
                            parts = line.split("|")
                            key = parts[1]
                            value = parts[2].split()[0]
                            data_dict[key] = value
                            data_list.append({"NameUn": key, "VerbalName": value})
        df_names = pd.DataFrame(data_list).drop_duplicates()
        protein_list = set(df_prot_eval_go.keys()).difference(set(df_names["NameUn"]))
        return df_names, protein_list

    @staticmethod
    def query_uniprot_for_names(
        protein_list: Set[str], base_url: str, df_names: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Queries UniProt for missing protein verbal names.
        """
        data_list = []
        for protein in protein_list:
            url = f"{base_url}{protein}.json"
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                data_list.append(
                    {
                        "primaryAccession": data.get("primaryAccession", ""),
                        "uniProtkbId": data.get("uniProtkbId", ""),
                    }
                )
        result_df = pd.DataFrame(data_list) if data_list else pd.DataFrame()
        result_df.rename(
            columns={"primaryAccession": "NameUn", "uniProtkbId": "VerbalName"},
            inplace=True,
        )
        return pd.concat([result_df, df_names], axis=0)

    @staticmethod
    def enrich_with_verbal_names(
        df_names: pd.DataFrame, df_ProtEvalGO: Dict[str, pd.DataFrame]
    ) -> Tuple[
        Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]
    ]:
        """
        Adds verbal names and separates GO vs. Ass dicts.
        """
        df_ProtEvalGO_upd, df_ProtEvalGO_updGO, df_ProtEvalGO_updAss = {}, {}, {}
        for key, df in df_ProtEvalGO.items():
            merged = pd.merge(
                df_names, df, left_on=["NameUn"], right_on=["Name"], how="right"
            )
            merged.drop(["NameUn"], axis=1, inplace=True)
            merged.rename(
                columns={"Name": "Name_Query", "VerbalName": "VerbalName_Query"},
                inplace=True,
            )
            df_ProtEvalGO_upd[key] = merged
            df_GO = merged[merged["BLASTvsAss_Flag"] != False]
            df_Ass = merged[merged["BLASTvsAss_Flag"] == False]
            if not df_GO.empty:
                df_ProtEvalGO_updGO[key] = df_GO
            if not df_Ass.empty:
                df_ProtEvalGO_updAss[key] = df_Ass
        return df_ProtEvalGO_upd, df_ProtEvalGO_updGO, df_ProtEvalGO_updAss

    @staticmethod
    def merge_with_pfp(
        input_dict: Dict[str, pd.DataFrame], df_main_100: Dict[str, pd.DataFrame]
    ) -> Dict[str, pd.DataFrame]:
        """
        Merges input dicts with PFP to enrich with PFP info.
        """
        dict_to_return = {}
        for prot in input_dict.keys():
            input_dict[prot].rename(columns={"Name_Query": "Name"}, inplace=True)
        for key, df in df_main_100.items():
            if key in input_dict.keys():
                merged = pd.merge(
                    df,
                    input_dict[key],
                    left_on=["Name", "GO_Term"],
                    right_on=["Name", "found_or_ass_go"],
                    how="inner",
                )
                merged.drop(
                    ["found_or_ass_go", "Actual_Blasted_GO_Term"], axis=1, inplace=True
                )
                if not merged.empty:
                    dict_to_return[key] = merged
        return dict_to_return

    @staticmethod
    def load_child_parent_relations(
        input_dir: str, df_main_100: Dict[str, pd.DataFrame]
    ) -> Dict[str, pd.DataFrame]:
        """
        Reads child-parent GO term relations.
        """
        df_ChPar = {}
        for filename in os.listdir(input_dir):
            name = filename[:-11]
            if name in df_main_100.keys() and filename.endswith(".csv"):
                df = pd.read_csv(
                    os.path.join(input_dir, filename),
                    dtype=str,
                    sep=",",
                    names=[
                        "Label",
                        "GO_TermChild",
                        "GO_TermParent",
                        "OverallParentalScore",
                    ],
                )
                df.drop(["Label", "OverallParentalScore"], axis=1, inplace=True)
                df["Name"] = name
                df_ChPar[name] = df
        return df_ChPar

    @staticmethod
    def enrich_with_parental_relations(
        dict_ChPar: Dict[str, pd.DataFrame],
        df_input: Dict[str, pd.DataFrame],
        df_main_100: Dict[str, pd.DataFrame],
    ) -> Dict[str, pd.DataFrame]:
        """
        Enriches with child-parent relationships and merges with PFP data.
        """
        dict_tmp, return_dict = {}, {}
        for key, df in dict_ChPar.items():
            if key in df_input.keys():
                merged = pd.merge(
                    df,
                    df_input[key],
                    left_on=["Name", "GO_TermChild"],
                    right_on=["Name", "Actual_Blasted_GO_Term"],
                    how="inner",
                )
                merged.drop(
                    ["Actual_Blasted_GO_Term", "found_or_ass_go"], axis=1, inplace=True
                )
                if not merged.empty:
                    dict_tmp[key] = merged
        for key, df in df_main_100.items():
            if key in dict_tmp.keys():
                merged = pd.merge(
                    df,
                    dict_tmp[key],
                    left_on=["Name", "GO_Term"],
                    right_on=["Name", "GO_TermParent"],
                    how="inner",
                )
                if not merged.empty:
                    return_dict[key] = merged
        return return_dict

    @staticmethod
    def organize_child_parent_percentages(file_path: str, output_dir: str) -> None:
        """
        Organizes child-parent percentages into a dictionary and saves it.
        """
        df_child_parent = {}
        df = pd.read_csv(
            file_path,
            sep=r"\s+",
            dtype=str,
            names=[
                "GO_Term_Child",
                "Parent",
                "Percentage",
                "a1",
                "a2",
                "a3",
                "a4",
                "a5",
            ],
        )
        df = df[["GO_Term_Child", "Parent", "Percentage"]]
        for term in set(df["GO_Term_Child"]):
            df_child_parent[term] = df[df["GO_Term_Child"] == term]
        saving_lists_of_locats_and_dicts([output_dir], [df_child_parent])

    @staticmethod
    def filter_by_self(
        input_dict: Dict[str, pd.DataFrame], columns_to_remove: List[str]
    ) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
        """
        Filters rows where query protein equals found protein.
        """
        filtered_dict, retained_dict = {}, {}
        for key, df in input_dict.items():
            filtered = df[df["VerbalName_Query"] != df["found_protein"]]
            retained = df[df["VerbalName_Query"] == df["found_protein"]]
            if not retained.empty:
                retained_dict[key] = retained
            filtered.drop(columns_to_remove, axis=1, inplace=True)
            filtered.drop_duplicates(inplace=True)
            if not filtered.empty:
                filtered_dict[key] = filtered
        return filtered_dict, retained_dict

    @staticmethod
    def apply_neglog10_transform(
        input_dict: Dict[str, pd.DataFrame]
    ) -> Dict[str, pd.DataFrame]:
        """
        Applies -log10(evalue) transformation.
        """
        edited_dict = {}
        for key, df in input_dict.items():
            df = df.copy()
            df["Num_found_evalue"] = df["found_evalue"].astype(float)
            if "GO_TermParent" in df.columns:
                df["Num_Percentage"] = df["Percentage"].astype(float)
                df["prelog+b_P"] = 2.097 - np.log10(df["Num_found_evalue"])
                df["-log+b_P"] = df["prelog+b_P"] * df["Num_Percentage"]
            else:
                df["-log+b"] = 2.097 - np.log10(df["Num_found_evalue"])
            edited_dict[key] = df.reset_index(drop=True)
        return edited_dict

    @staticmethod
    def validate_single_flag(
        input_dict: Dict[str, pd.DataFrame], column_name: str, expected_value: Any
    ) -> None:
        """
        Validates whether the unique value in the column matches the expected value.
        """
        df = pd.concat(input_dict.values(), axis=0, ignore_index=True)
        if column_name not in df.columns:
            raise KeyError(f"Column '{column_name}' not found.")
        unique_values = df[column_name].unique()
        if len(unique_values) != 1:
            raise ValueError(f"Multiple values found: {unique_values}")
        actual = unique_values[0]
        if bool(actual) != bool(expected_value):
            raise ValueError(f"Expected {expected_value} but found {actual}.")
        print(
            f"[CHECKPOINT PASSED for '{column_name}'] - Found {actual}, matches expected {expected_value}"
        )

    @staticmethod
    def checkpoint_results(
        df_main_100: Dict[str, pd.DataFrame],
        df_main_allCP_GO: Dict[str, pd.DataFrame],
        df_main_allCP_Ass: Dict[str, pd.DataFrame],
        df_main_allGO: Dict[str, pd.DataFrame],
        df_main_allAss: Dict[str, pd.DataFrame],
    ) -> None:
        """
        Checks whether all original results are retrieved after merging.
        """
        columns = [
            "GO_Term",
            "Ontology",
            "Out_of_Total",
            "Raw_Score",
            "Rank",
            "Verbal_Desc",
            "Name",
        ]

        merged = pd.concat(
            [
                pd.concat(df.values(), axis=0)[columns]
                for df in [
                    df_main_allCP_GO,
                    df_main_allCP_Ass,
                    df_main_allGO,
                    df_main_allAss,
                ]
            ],
            axis=0,
        ).drop_duplicates()

        original = pd.concat(df_main_100.values(), axis=0)

        compare = original.merge(merged, on=columns, how="outer", indicator=True)
        both = compare[compare["_merge"] == "both"]

        if (
            both[columns]
            .sort_values(columns)
            .reset_index(drop=True)
            .equals(original[columns].sort_values(columns).reset_index(drop=True))
        ):
            print("[CHECKPOINT PASSED] - All original results retrieved and match.")
        else:
            missing = original.shape[0] - both.shape[0]
            raise ValueError(
                f"[CHECKPOINT FAILED] - {missing} rows missing after merging."
            )

        del merged, original, compare, both
        gc.collect()

    @staticmethod
    def fast_merge_per(
        input_dict: Dict[str, pd.DataFrame],
        child_parent_dict: Dict[str, pd.DataFrame],
    ) -> Dict[str, pd.DataFrame]:
        """
        Efficiently merges child-parent relationships across proteins.
        """
        result_dict = {}
        for key, df in input_dict.items():
            child_terms = set(df["GO_TermChild"])
            tmp_dfs = [
                child_parent_dict[term]
                for term in child_terms
                if term in child_parent_dict
            ]
            tmp_full = pd.concat(tmp_dfs, axis=0, ignore_index=True)
            merged = pd.merge(
                df,
                tmp_full,
                left_on=["GO_TermChild", "GO_TermParent"],
                right_on=["GO_Term_Child", "Parent"],
                how="inner",
            )
            if not merged.empty:
                result_dict[key] = merged
        return result_dict
