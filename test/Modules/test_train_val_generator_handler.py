import tempfile
import unittest
from unittest.mock import patch, MagicMock

import numpy as np
import pandas as pd

from Modules.train_val_generator_handler import TrainValidationSplitHandler


class TestTrainValidationSplitHandler(unittest.TestCase):

    @patch("os.listdir")
    @patch("pandas.read_parquet")
    def test_load_4_class_dictionaries(self, mock_read_parquet, mock_listdir):
        mock_listdir.side_effect = [
            ["a.parquet"],
            ["b.parquet"],
            ["c.parquet"],
            ["d.parquet"],
        ]
        mock_read_parquet.return_value = pd.DataFrame({"col": [1]})

        result = TrainValidationSplitHandler.load_4_class_dictionaries(
            "a", "b", "c", "d"
        )

        self.assertEqual(len(result), 4)
        for subdict in result.values():
            self.assertTrue(isinstance(list(subdict.values())[0], pd.DataFrame))

    def test_combine_4_class_dictionaries(self):
        df_a = pd.DataFrame({"Name": ["A1"], "-log+b": [0.1]})
        df_b = pd.DataFrame({"Name": ["B1"], "-log+b": [0.2]})
        df_c = pd.DataFrame(
            {
                "Name": ["C1"],
                "GO_TermChild": ["x"],
                "GO_TermParent": ["y"],
                "Percentage": [0.5],
                "Num_Percentage": [0.1],
                "prelog+b_P": [0.3],
                "-log+b_P": [0.4],
            }
        )
        df_d = pd.DataFrame(
            {
                "Name": ["D1"],
                "GO_TermChild": ["x"],
                "GO_TermParent": ["y"],
                "Percentage": [0.5],
                "Num_Percentage": [0.1],
                "prelog+b_P": [0.3],
                "-log+b_P": [0.4],
            }
        )

        dict_a = {"A": df_a}
        dict_b = {"B": df_b}
        dict_c = {"C": df_c}
        dict_d = {"D": df_d}

        combined = TrainValidationSplitHandler.combine_4_class_dictionaries(
            dict_a, dict_b, dict_c, dict_d
        )
        self.assertEqual(combined.shape[0], 4)
        self.assertIn("Class", combined.columns)

    @patch("pandas.DataFrame.to_parquet")
    def test_perform_group_shuffle_split(self, mock_to_parquet):
        df = pd.DataFrame(
            {"Name": [f"P{i}" for i in range(10)] * 5, "Feature": np.random.rand(50)}
        )

        log_mock = MagicMock()

        with tempfile.TemporaryDirectory() as tmpdirname:
            TrainValidationSplitHandler.perform_group_shuffle_split(
                log_handler=log_mock,
                combined_df=df,
                output_dir=tmpdirname,
                n_splits_outer=2,
                n_splits_inner=1,
            )

            # Check that parquet files were saved
            self.assertTrue(mock_to_parquet.called)
            self.assertGreaterEqual(
                mock_to_parquet.call_count, 6
            )  # 2 outer * 3 inner (test, train_val, train, val)


if __name__ == "__main__":
    unittest.main()
