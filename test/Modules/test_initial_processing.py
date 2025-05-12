import unittest
from unittest.mock import patch, mock_open, MagicMock

import pandas as pd

from Modules.initial_preprocessing import InitialPreprocessor


class TestInitialPreprocessor(unittest.TestCase):

    @patch("os.listdir")
    @patch("pandas.read_csv")
    def test_load_pfp_predictions(self, mock_read_csv, mock_listdir):
        mock_listdir.return_value = ["protein1.csv"]
        df = pd.DataFrame(
            {
                "GO_Term": ["GO:0001"],
                "Ontology": ["BP"],
                "Out_of_Total": ["10/100"],
                "Raw_Score": [150.0],
                "Rank": ["1"],
                "Verbal_Desc": ["desc"],
            }
        )
        mock_read_csv.return_value = df
        result = InitialPreprocessor.load_pfp_predictions("mock_dir")
        self.assertIn("protein1", result)

    @patch("glob.glob")
    def test_extract_verbal_names(self, mock_glob):
        mock_glob.side_effect = [["/mock_dir"], ["/mock_dir/blast/protein_blast"]]
        file_content = "sp|P12345|TestProtein description\n"
        with patch("builtins.open", mock_open(read_data=file_content)):
            df_ProtEvalGO = {"P12345": pd.DataFrame()}
            df_names, protein_list = InitialPreprocessor.extract_verbal_names(
                "~", df_ProtEvalGO
            )
            self.assertIn("P12345", df_names["NameUn"].values)

    @patch("requests.get")
    def test_query_uniprot_for_names(self, mock_get):
        response = MagicMock()
        response.status_code = 200
        response.json.return_value = {
            "primaryAccession": "P12345",
            "uniProtkbId": "TestProtein",
        }
        mock_get.return_value = response
        df_names = pd.DataFrame()
        result = InitialPreprocessor.query_uniprot_for_names(
            {"P12345"}, "https://mock.url/", df_names
        )
        self.assertIn("P12345", result["NameUn"].values)

    def test_enrich_with_verbal_names(self):
        df_names = pd.DataFrame(
            {"NameUn": ["protein1"], "VerbalName": ["Protein Name"]}
        )
        df_prot_eval_go = {
            "protein1": pd.DataFrame(
                {
                    "Name": ["protein1"],
                    "Type_of_GO": ["TypeisGO"],
                    "BLASTvsAss_Flag": [True],
                    "found_or_ass_go": ["GO:0001"],
                    "Actual_Blasted_GO_Term": ["GO:0001"],
                    "found_protein": ["proteinA"],
                    "found_evalue": ["1e-5"],
                }
            )
        }

        upd, updGO, updAss = InitialPreprocessor.enrich_with_verbal_names(
            df_names, df_prot_eval_go
        )

        # Basic correctness checks
        self.assertIn("protein1", upd)
        self.assertIn("protein1", updGO)
        # Since BLASTvsAss_Flag is True, Ass should be empty
        self.assertNotIn("protein1", updAss)

    def test_merge_with_pfp(self):
        df_input = {
            "protein1": pd.DataFrame(
                {
                    "Name_Query": ["protein1"],
                    "found_or_ass_go": ["GO:0001"],
                    "Actual_Blasted_GO_Term": ["GO:0001"],
                }
            )
        }
        df_main_100 = {
            "protein1": pd.DataFrame({"Name": ["protein1"], "GO_Term": ["GO:0001"]})
        }
        result = InitialPreprocessor.merge_with_pfp(df_input, df_main_100)
        self.assertIn("protein1", result)

    def test_filter_by_self(self):
        df = pd.DataFrame(
            {
                "VerbalName_Query": ["A", "B"],
                "found_protein": ["A", "C"],
                "extra": [1, 2],
            }
        )
        input_dict = {"protein": df}
        filtered, retained = InitialPreprocessor.filter_by_self(input_dict, ["extra"])
        self.assertIn("protein", filtered)
        self.assertIn("protein", retained)

    def test_apply_neglog10_transform(self):
        df = pd.DataFrame({"found_evalue": [1e-5]})
        input_dict = {"protein": df}
        result = InitialPreprocessor.apply_neglog10_transform(input_dict)
        self.assertIn("-log+b", result["protein"].columns)

    def test_validate_single_flag(self):
        input_dict = {"protein": pd.DataFrame({"flag": [True]})}
        InitialPreprocessor.validate_single_flag(input_dict, "flag", True)

    def test_fast_merge_per(self):
        input_dict = {
            "protein": pd.DataFrame(
                {"GO_TermChild": ["GO:0002"], "GO_TermParent": ["GO:0001"]}
            )
        }
        child_parent_dict = {
            "GO:0002": pd.DataFrame(
                {"GO_Term_Child": ["GO:0002"], "Parent": ["GO:0001"]}
            )
        }
        result = InitialPreprocessor.fast_merge_per(input_dict, child_parent_dict)
        self.assertIn("protein", result)


if __name__ == "__main__":
    unittest.main()
