import glob
import json
import random
from collections import defaultdict
from itertools import chain
from typing import Dict, Tuple

import networkx as nx
import numpy as np

from Utils.utils import *
from pygosemsim import similarity


class PyGOsemsimHandler:
    """
    A handler class for all pygosemsim-based data preparation and analysis utilities.
    All original procedural functions have been encapsulated here as static methods.
    """

    @staticmethod
    def get_dataframe(
        single_predict_df: pd.DataFrame,
    ) -> Tuple[pd.DataFrame, List[str]]:
        single_predict_df = single_predict_df.reset_index(drop=True)
        single_predict_prot_list = list(single_predict_df["Name"].drop_duplicates())
        single_predict_df["Ontology"] = single_predict_df["Ontology"].str.upper()
        single_predict_df["PFP_Pred_GO_term"] = "GO:" + single_predict_df["GO"].astype(
            str
        )
        single_predict_df = single_predict_df.drop(["GO"], axis=1)
        single_predict_df.rename(
            columns={"Name": "Protein", "Ontology": "Predicted_Ontology"}, inplace=True
        )
        print("length of the protein set", len(single_predict_prot_list))
        return single_predict_df, single_predict_prot_list

    @staticmethod
    def get_6_species_with_gaf(
        updated_uniprot_6_species_with_g_oterms_gaf: pd.DataFrame,
    ) -> pd.DataFrame:
        updated_uniprot_6_species_with_g_oterms_gaf_deduped = (
            updated_uniprot_6_species_with_g_oterms_gaf.drop_duplicates()
        )
        updated_uniprot_6_species_with_g_oterms_gaf1 = (
            updated_uniprot_6_species_with_g_oterms_gaf_deduped[
                ["Protein", "GO_term", "Ontology"]
            ]
        )
        updated_uniprot_6_species_with_g_oterms_gaf2 = (
            updated_uniprot_6_species_with_g_oterms_gaf1.copy()
        )
        updated_uniprot_6_species_with_g_oterms_gaf2.rename(
            columns={"GO_term": "GAF_GO_term"}, inplace=True
        )
        return updated_uniprot_6_species_with_g_oterms_gaf2

    @staticmethod
    def combining_partial_files(
        file_path: str, sw_path: str
    ) -> Tuple[defaultdict, List, List]:
        list_gaf = glob.glob(file_path + "/annot_dic*")
        file_path = os.path.abspath(file_path)

        count = 1
        list_var_partial = []
        for item in list_gaf:
            with open(item, "r") as f:
                var = json.load(f)
                count += 1
                list_var_partial.append(var)

        with open(file_path + "/y_combined.json", "r") as f:
            y_comb = json.load(f)

        simple_dic = defaultdict(set)
        for dic in list_var_partial:
            for protein, details in dic.items():
                simple_dic[protein].update(details["annotation"].keys())

        with open(sw_path, "r") as f:
            swiss_prot_set = set(line.split("\t")[0] for line in f)

        y_comb1 = [y for y in y_comb if y is not None and y[0] in swiss_prot_set]
        simple_y_comb = defaultdict(list)
        for item in y_comb1:
            if (
                item[2] not in {"contributes_to", "colocalizes_with"}
                and "NOT" not in item[2]
            ):
                simple_y_comb[item[0]] = item[1]

        for key in simple_y_comb:
            simple_dic[key].add(simple_y_comb[key])

        keys_to_delete = [key for key in simple_dic if simple_dic[key] == set()]
        for key in keys_to_delete:
            del simple_dic[key]

        go_terms = []
        proteins = []
        for protein in simple_dic:
            proteins.append(protein)
            for go in simple_dic[protein]:
                go_terms.append(go)

        return simple_dic, go_terms, proteins

    @staticmethod
    def propagate_until_root(
        G: nx.DiGraph,
        alt_dic: dict,
        repl_dic: dict,
        deleted_list: list,
        gaf_df: pd.DataFrame,
    ) -> pd.DataFrame:
        gaf_df["GAF_GO_term"] = gaf_df["GAF_GO_term"].replace(alt_dic)
        gaf_df["GAF_GO_term"] = gaf_df["GAF_GO_term"].replace(repl_dic)
        gaf_df = gaf_df[~gaf_df["GAF_GO_term"].isin(deleted_list)]

        grouped = gaf_df.groupby("Protein")["GAF_GO_term"].apply(list)
        expanded_annotations = {}
        for prot, terms in grouped.items():
            new_terms = set()
            for go in terms:
                if go in G:
                    ancestors = nx.ancestors(G, go)
                    new_terms.update(ancestors)
                    new_terms.add(go)
            expanded_annotations[prot] = list(new_terms)

        exploded = [
            (protein, go) for protein, gos in expanded_annotations.items() for go in gos
        ]
        return pd.DataFrame(exploded, columns=["Protein", "GAF_GO_term"])

    @staticmethod
    def sort_real_go_dict(real_go_dic: Dict[str, List[str]]) -> Dict[str, List[str]]:
        return {prot: sorted(go_list) for prot, go_list in real_go_dic.items()}

    @staticmethod
    def time_diff_check_1(
        association_df: pd.DataFrame,
        category_df: pd.DataFrame,
        child_parent_percentages_df: pd.DataFrame,
        go_term_1: pd.DataFrame,
        go_term_2: pd.DataFrame,
    ) -> List[str]:
        set1 = set(association_df["GO1"].to_list() + association_df["GO2"].to_list())
        set2 = set(category_df["GO"].to_list())
        set3 = set(
            child_parent_percentages_df["Child"].to_list()
            + child_parent_percentages_df["Parent"].to_list()
        )
        set4 = set(go_term_1["Id"].to_list() + go_term_2["Id"].to_list())
        return list(set1.union(set2).union(set3).union(set4))

    @staticmethod
    def time_diff_check_2(
        G: nx.DiGraph,
        loc_2019_basic: str,
        list_combined_pfp: List[str],
        repl_dic: dict,
        alt_dic: dict,
    ) -> List[str]:
        g_nodes = set(G.nodes())
        with open(loc_2019_basic, "r") as f:
            go_terms_2019 = set(line.strip() for line in f.readlines())
        not_found = go_terms_2019 - set(list_combined_pfp)
        not_found = not_found - set(repl_dic.keys())
        not_found = not_found - set(alt_dic.keys())
        not_found = [go for go in not_found if go in g_nodes]
        return not_found

    @staticmethod
    def time_diff_check_3_update_real_go_dict1_part1(
        realGO_dic_sorted: Dict[str, List[str]],
        del_dict: Dict[str, List[str]],
        add_dict: Dict[str, List[str]],
        not_found: List[str],
    ) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[str, List[str]]]:
        for prot, go_list in realGO_dic_sorted.items():
            to_remove = [go for go in go_list if go in not_found]
            if to_remove:
                del_dict[prot] = to_remove
            add_dict[prot] = []
        return realGO_dic_sorted, del_dict, add_dict

    @staticmethod
    def time_diff_check_3_update_real_go_dict1_part2(
        realGO_dic_sorted: Dict[str, List[str]], not_found: List[str]
    ) -> Tuple[Dict[str, List[str]], List[str]]:
        list_removed = []
        for prot, go_list in realGO_dic_sorted.items():
            new_list = [go for go in go_list if go not in not_found]
            if len(new_list) != len(go_list):
                list_removed.extend(set(go_list) - set(new_list))
            realGO_dic_sorted[prot] = new_list
        return realGO_dic_sorted, list_removed

    @staticmethod
    def propagate_gaf_until_root(
        simple_dic_var: defaultdict,
        G: nx.DiGraph,
        alt_dic: dict,
        repl_dic: dict,
        deleted_list: list,
    ) -> Dict[str, List[str]]:
        propagated_dic = {}
        for prot, annots in simple_dic_var.items():
            go_terms = annots["annotation"].keys()
            expanded_terms = set()
            for go in go_terms:
                go = repl_dic.get(alt_dic.get(go, go), go)
                if go not in deleted_list and go in G:
                    ancestors = nx.ancestors(G, go)
                    expanded_terms.update(ancestors)
                    expanded_terms.add(go)
            propagated_dic[prot] = list(expanded_terms)
        return propagated_dic

    @staticmethod
    def breaking_propagation_sw_to_ont(
        propagation_dic: Dict[str, List[str]]
    ) -> Tuple[List[List[str]], List[List[str]], List[List[str]]]:
        BP, CC, MF = [], [], []
        for prot, terms in propagation_dic.items():
            for go in terms:
                if go.startswith("GO:0008"):
                    BP.append([prot, go])
                elif go.startswith("GO:0005"):
                    CC.append([prot, go])
                elif go.startswith("GO:0003"):
                    MF.append([prot, go])
        return BP, CC, MF

    @staticmethod
    def joined_lsit_fo_function(
        joined_list: List[List[str]], label: str
    ) -> pd.DataFrame:
        return pd.DataFrame(joined_list, columns=["Protein", "GO"]).assign(Label=label)

    @staticmethod
    def dic_with_count(G: nx.DiGraph, df: pd.DataFrame) -> Dict[str, int]:
        counts = defaultdict(int)
        for go in df["GO"]:
            counts[go] += 1
        for node in G.nodes():
            counts[node] = counts.get(node, 0)
        return counts

    @staticmethod
    def get_ic(count_dic: Dict[str, int], root: str) -> Dict[str, float]:
        total = count_dic[root]
        return {
            term: -np.log(count / total) if count > 0 else 0.0
            for term, count in count_dic.items()
        }

    @staticmethod
    def checkpoint1(
        dic_bp: Dict[str, float], dic_cc: Dict[str, float], dic_mf: Dict[str, float]
    ) -> int:
        set_bp = set(dic_bp.keys())
        set_cc = set(dic_cc.keys())
        set_mf = set(dic_mf.keys())
        return len(set_bp & set_cc) + len(set_bp & set_mf) + len(set_cc & set_mf)

    def edit_and_shape_by_ontology_pred(
        self, input_list: list, alt_dict: dict, repl_dict: dict, deleted_list: list
    ) -> list:
        """
        Edits GO terms predicted list using alternative, replacement and deleted lists.
        """
        edited_list_alt = [alt_dict.get(go, go) for go in input_list]
        edited_list_repl = [repl_dict.get(go, go) for go in edited_list_alt]
        edited_list_deleted = [go for go in edited_list_repl if go not in deleted_list]
        return edited_list_deleted

    def edit_and_shape_by_ontology_real(
        self, input_list: list, alt_dict: dict, repl_dict: dict, deleted_list: list
    ) -> list:
        """
        Edits GO terms real list using alternative, replacement and deleted lists.
        """
        combined = []
        for sublist in input_list:
            edited_list_alt = [alt_dict.get(go, go) for go in sublist]
            edited_list_repl = [repl_dict.get(go, go) for go in edited_list_alt]
            edited_list_deleted = [
                go for go in edited_list_repl if go not in deleted_list
            ]
            combined.append(edited_list_deleted)
        if input_list != combined:
            print("False")
        return combined

    @staticmethod
    def edit_and_prepare_binary(
        input_list: list, alt_dict: dict, repl_dict: dict, deleted_list: list
    ) -> list:
        """
        Edits and prepares GO terms list for binary scoring.
        """
        edited_list_alt = [alt_dict.get(go, go) for go in input_list]
        edited_list_repl = [repl_dict.get(go, go) for go in edited_list_alt]
        edited_list_deleted = [go for go in edited_list_repl if go not in deleted_list]
        return edited_list_deleted

    @staticmethod
    def create_binary(predicted_list: list, real_list: list) -> tuple:
        """
        Generates binary labels between predicted and real GO terms.
        """
        binary = []
        real_missing = []
        predicted_set = list(set(predicted_list))

        for pred in predicted_set:
            if pred in real_list:
                binary.append([pred, pred, 1])
            else:
                binary.append([pred, "Not Found", 0])

        for real in real_list:
            if real not in predicted_set:
                real_missing.append([real, "Real Missing from Pred"])

        return binary, real_missing

    @staticmethod
    def make_lin_tuples(predicted_list: list, real_nested_list: list) -> list:
        """
        Makes Lin tuples from predicted and real GO term lists.
        """
        from itertools import chain

        list_lin_tuples = []
        flattened_real = list(set(list(chain.from_iterable(real_nested_list))))

        initial_elements = [sublist[0] for sublist in real_nested_list if sublist]
        propagated_elements = [sublist[1:] for sublist in real_nested_list if sublist]

        for pred in set(predicted_list):
            if not real_nested_list:
                list_lin_tuples.append(
                    [
                        pred,
                        "Placeholder",
                        "No ground truth terms/No term in same ontology",
                    ]
                )
                continue

            if pred in initial_elements:
                list_lin_tuples.append([pred, pred, "node", "found"])
            elif any(pred in sub for sub in propagated_elements):
                list_lin_tuples.append([pred, pred, "Parent node", "found"])
            else:
                for real in flattened_real:
                    list_lin_tuples.append(
                        [pred, real, "general pool of nodes", "not found"]
                    )

        return list_lin_tuples

    @staticmethod
    def generate_binary_dataframe(protein_ontologies: dict) -> pd.DataFrame:
        """
        Generates DataFrame of binary scores from protein ontologies.
        """
        from collections import defaultdict

        binary_dict = defaultdict(list)
        for protein in protein_ontologies:
            for pair in protein_ontologies[protein]["binary"]:
                binary_dict[protein].append([pair[0], pair[2]])
        dic1 = {
            key: pd.DataFrame(value, columns=["PredictedGO", "BinaryScore"]).assign(
                Protein=key
            )
            for key, value in binary_dict.items()
        }
        binary_df = pd.concat(
            dic1.values(), axis=0, ignore_index=True
        ).drop_duplicates()
        return binary_df

    @staticmethod
    def cleaning_predictions(
        df: pd.DataFrame, alt_dict: dict, repl_dict: dict, deleted_list: list
    ) -> pd.DataFrame:
        """
        Cleans the predicted dataframe using mappings for alt/replacement/deleted terms.
        """

        df["changed to altD"] = df["PFP_Pred_GO_term"].apply(
            lambda x: alt_dict.get(x, x)
        )
        df["presence of altD"] = df["PFP_Pred_GO_term"].apply(
            lambda x: int(x in alt_dict)
        )
        df["changed to replD"] = df["changed to altD"].apply(
            lambda x: repl_dict.get(x, x)
        )
        df["presence of replD"] = df["changed to altD"].apply(
            lambda x: int(x in repl_dict)
        )

        deleted_dict = dict.fromkeys(deleted_list, None)
        df["changed to deleted"] = df["changed to replD"].apply(
            lambda x: "Deleted" if x in deleted_dict else x
        )
        df["presence of deleted"] = df["changed to replD"].apply(
            lambda x: int(x in deleted_dict)
        )

        return df[df["presence of deleted"] != 1].copy(deep=True)

    @staticmethod
    def query_uniprot_for_stats(base_url: str, protein_list: list) -> pd.DataFrame:
        """
        Queries Uniprot REST API for missing protein features.
        """
        import requests

        data_list = []
        for protein in protein_list:
            url = f"{base_url}{protein}.json"
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                data_list.append(
                    {
                        "Name": data.get("primaryAccession", ""),
                        "AnnotationScore": data.get("annotationScore", ""),
                        "ProteinExistence": data.get("proteinExistence", "")[0],
                        "Seq_length": data.get("sequence", {}).get("length", ""),
                    }
                )
            else:
                print(f"Error fetching data for {protein}: {response.status_code}")

        return pd.DataFrame(data_list) if data_list else pd.DataFrame()

    @staticmethod
    def add_go_term_depth(graph_obj, df: pd.DataFrame) -> pd.DataFrame:
        """
        Adds GO term depth to dataframe based on GO graph structure.
        """
        import networkx as nx

        go_list = list(set(df["Updated_PFP_Pred_GO_term"].to_list()))
        ontology_dict = dict(
            zip(df["Updated_PFP_Pred_GO_term"], df["Predicted_Ontology_Retained"])
        )
        depth_records = []
        for go in go_list:
            if go in graph_obj.nodes():
                root = {"C": "GO:0005575", "F": "GO:0003674", "P": "GO:0008150"}[
                    ontology_dict[go]
                ]
                depth = nx.shortest_path_length(graph_obj, source=root, target=go)
                depth_records.append({"Updated_PFP_Pred_GO_term": go, "Depth": depth})

        depth_df = pd.DataFrame(depth_records)
        df_merged = pd.merge(df, depth_df, on="Updated_PFP_Pred_GO_term", how="inner")
        df_renamed = df_merged.rename(
            columns={
                "filtered_count_AA": "DirectBlastAA",
                "filtered_count_BB": "AssociatedBlastBB",
                "filtered_count_CC": "DirectParentCC",
                "filtered_count_DD": "AssociatedParentDD",
            }
        )

        exclude_terms = ["GO:0005575", "GO:0003674", "GO:0008150"]
        return df_renamed[~df_renamed["Original_PFP_Pred_GO_term"].isin(exclude_terms)]

    def get_ont_lists(
        self,
        protein: str,
        real_go_dict: dict,
        pred_go_dict: dict,
        alt_dict: dict,
        repl_dict: dict,
        deleted_list: list,
    ) -> tuple:
        """
        Generates lists of tuples per ontology and corresponding missing dataset.
        """
        temp_r_c, temp_r_f, temp_r_p = [], [], []
        temp_p_c, temp_p_f, temp_p_p = [], [], []
        temp_p2 = []
        temp_r1 = real_go_dict[protein]
        temp_r2 = [temp_r1[x][1] for x in temp_r1]
        temp_r2_flat = list(set(list(chain.from_iterable(temp_r2))))

        if protein in real_go_dict:
            temp_r1 = real_go_dict[protein]
            temp_r_c, temp_r_f, temp_r_p = [], [], []
            for go, ont in temp_r1.items():
                if ont[0] == "C":
                    temp_r_c.append(ont[1])
                if ont[0] == "F":
                    temp_r_f.append(ont[1])
                if ont[0] == "P":
                    temp_r_p.append(ont[1])

        if protein in pred_go_dict:
            temp_p1 = pred_go_dict[protein]
            temp_p2 = [x[0] for x in temp_p1]
            for each in temp_p1:
                if each[1] == "C":
                    temp_p_c.append(each[0])
                if each[1] == "F":
                    temp_p_f.append(each[0])
                if each[1] == "P":
                    temp_p_p.append(each[0])

        edited_temp_p_c = self.edit_and_shape_by_ontology_pred(
            temp_p_c, alt_dict, repl_dict, deleted_list
        )
        edited_temp_p_f = self.edit_and_shape_by_ontology_pred(
            temp_p_f, alt_dict, repl_dict, deleted_list
        )
        edited_temp_p_p = self.edit_and_shape_by_ontology_pred(
            temp_p_p, alt_dict, repl_dict, deleted_list
        )

        edited_temp_r_c = self.edit_and_shape_by_ontology_real(
            temp_r_c, alt_dict, repl_dict, deleted_list
        )
        edited_temp_r_f = self.edit_and_shape_by_ontology_real(
            temp_r_f, alt_dict, repl_dict, deleted_list
        )
        edited_temp_r_p = self.edit_and_shape_by_ontology_real(
            temp_r_p, alt_dict, repl_dict, deleted_list
        )

        edited_temp_p = self.edit_and_prepare_binary(
            temp_p2, alt_dict, repl_dict, deleted_list
        )
        edited_temp_r = self.edit_and_prepare_binary(
            temp_r2_flat, alt_dict, repl_dict, deleted_list
        )

        binary, real_missing = self.create_binary(edited_temp_p, edited_temp_r)

        list_lin_tuples_c = self.make_lin_tuples(edited_temp_p_c, edited_temp_r_c)
        list_lin_tuples_f = self.make_lin_tuples(edited_temp_p_f, edited_temp_r_f)
        list_lin_tuples_p = self.make_lin_tuples(edited_temp_p_p, edited_temp_r_p)

        return (
            edited_temp_p_c,
            edited_temp_r_c,
            binary,
            real_missing,
            list_lin_tuples_c,
            list_lin_tuples_f,
            list_lin_tuples_p,
            temp_r2_flat,
            temp_r_c,
            temp_r_f,
            temp_r_p,
            temp_p2,
            temp_p_c,
            temp_p_f,
            temp_p_p,
            protein,
        )

    @staticmethod
    def add_root_flag(
        input_df: pd.DataFrame, root_term: str, ontology_label: str
    ) -> pd.DataFrame:
        """
        Adds a root flag to indicate if GO term equals root term and assigns ontology label.
        """
        root_flags = [
            1 if go == root_term else 0
            for go in input_df["PFP_Predicted_and_edited"].to_list()
        ]
        input_df["Root_Flag"] = root_flags
        input_df["LinOnt1"] = ontology_label
        return input_df

    def choosing_highest_score_within_pred_real(
        self,
        graph_obj,
        ic_dict: dict,
        protein_ontologies: dict,
        ontology_label: str,
        not_found_list: list,
    ) -> dict:
        """
        Calculates highest Lin score between predicted and real GO terms for all proteins.
        """
        return_dict = defaultdict(list)
        for prot in protein_ontologies.keys():
            tmp_d, result, result1, one_lin_result, result2 = self.lin_score(
                graph_obj,
                ic_dict,
                prot,
                ontology_label,
                protein_ontologies,
                not_found_list,
            )
            return_dict[prot] = result2
        return return_dict

    @staticmethod
    def lin_score(
        graph_obj,
        ic_dict: dict,
        protein: str,
        ontology_label: str,
        protein_ontologies: dict,
        not_found_list: list,
    ) -> tuple:
        """
        Calculates Lin similarity score for a protein's GO terms.
        """
        lin_score_list = []
        similarity.precalc_lower_bounds(graph_obj)
        go_pair_list = protein_ontologies[protein][ontology_label]
        if [
            pair for pair in go_pair_list if pair[1] not in not_found_list
        ] != go_pair_list:
            go_pair_list = [
                pair for pair in go_pair_list if pair[1] not in not_found_list
            ]

        for pair in go_pair_list:
            if pair[0] == pair[1]:
                lin_score_list.append(
                    [1, pair[0], pair[1], "Both GO terms are the same"]
                )
                continue
            if pair[1] == "Placeholder":
                lin_score_list.append(
                    [
                        0,
                        pair[0],
                        pair[1],
                        "No term in same ontology in found in Ground Truth",
                    ]
                )
                continue

            mica = similarity.lowest_common_ancestor(graph_obj, pair[0], pair[1])
            if mica is None:
                lin_score_list.append(
                    ["MICA is None", pair[0], pair[1], "One of the GO terms not found"]
                )
                continue
            lin_score = 2 * ic_dict[mica] / (ic_dict[pair[0]] + ic_dict[pair[1]])
            lin_score_list.append([abs(lin_score), pair, "Normal"])

        result_temp = defaultdict(list)
        one_lin_result = defaultdict(list)
        for item in lin_score_list:
            if item[0] != 1 and item[0] != "MICA is None" and item[2] != "Placeholder":
                result_temp[item[1][0]].append([item[0], item[1][0], item[1][1]])
            if item[0] == 0 and item[2] == "Placeholder":
                one_lin_result[item[1]].append(
                    [0, item[1], "No terms in this ontology"]
                )
            if item[0] == 1:
                one_lin_result[item[1]].append([1, item[1], item[1]])

        result = {}
        for key, values in result_temp.items():
            max_val = max(values, key=lambda x: x[0])[0]
            result[key] = [
                (val[0], val[1], val[2]) for val in values if val[0] == max_val
            ]
        result1 = {}
        for key, values in result.items():
            if len(values) > 1:
                result1[key] = random.choice(values)
            else:
                result1[key] = values[0]
        result2 = result1.copy()
        result2.update(one_lin_result)
        return result_temp, result, result1, one_lin_result, result2

    @staticmethod
    def lin_dic_to_df(input_dic: dict) -> pd.DataFrame:
        """
        Converts Lin similarity results dictionary to DataFrame.
        """
        rows, rows1, rows_n = [], [], []
        for prot, values in input_dic.items():
            for v1, v2 in values.items():
                if len(v2) == 1 and v2[0][0] == 1:
                    rows1.append(
                        {
                            "LinScore": v2[0][0],
                            "PFP_Predicted_and_edited": v2[0][1],
                            "PFP_Related/NearestNeighbor": v2[0][1],
                            "Protein": prot,
                        }
                    )
                elif len(v2) == 1 and v2[0][0] == 0:
                    rows_n.append(
                        {
                            "LinScore": v2[0][0],
                            "PFP_Predicted_and_edited": v2[0][1],
                            "PFP_Related/NearestNeighbor": "No terms in this ontology",
                            "Protein": prot,
                        }
                    )
                else:
                    rows.append(
                        {
                            "LinScore": v2[0],
                            "PFP_Predicted_and_edited": v2[1],
                            "PFP_Related/NearestNeighbor": v2[2],
                            "Protein": prot,
                        }
                    )
        df = pd.DataFrame(rows + rows1 + rows_n)
        return df

    @staticmethod
    def reshape_real_and_pred_propagated_go_terms(
        single_predict_df: pd.DataFrame, uniprot_propagated_gaf2_df: pd.DataFrame
    ) -> tuple:
        """
        Generates real GO term dictionary from propagated GAF2 file and predicted GO dictionary from single_predic_df.
        """
        pred_go_dict = defaultdict(list)
        for index, row in single_predict_df.iterrows():
            protein = row["Protein"]
            pred_go_dict[protein].append(
                (row["PFP_Pred_GO_term"], row["Predicted_Ontology"])
            )

        real_go_dict = {}
        for index, row in uniprot_propagated_gaf2_df.iterrows():
            protein = row["Protein"]
            go_term = row["GAF_GO_term"]
            propagated_terms = row["propagated GO Terms"]
            ontology = row["Ontology"]

            if protein not in real_go_dict:
                real_go_dict[protein] = {}
            real_go_dict[protein][go_term] = [ontology, propagated_terms]

        return real_go_dict, pred_go_dict
