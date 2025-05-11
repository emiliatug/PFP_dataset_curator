import logging
import os
from typing import Dict

import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

logger = logging.getLogger(__name__)


class TrainValidationSplitHandler:
    """
    A class that handles reading, combining, and splitting of dataset dictionaries
    for train-validation-test splits using GroupShuffleSplit.
    """

    def __init__(self):
        pass

    @staticmethod
    def load_4_class_dictionaries(
        dir_class_a: str, dir_class_b: str, dir_class_c: str, dir_class_d: str
    ) -> Dict[str, Dict[str, pd.DataFrame]]:
        """
        Reads 4 directories containing parquet files and creates 4 dictionaries.

        Returns:
            Dictionary containing 4 sub-dictionaries of DataFrames.
        """
        directories = {
            "df_main_allGO_hist0": dir_class_a,
            "df_main_allAss_hist0": dir_class_b,
            "df_main_allCP_GO_hist0": dir_class_c,
            "df_main_allCP_Ass_hist0": dir_class_d,
        }

        df_class_a, df_class_b, df_class_c, df_class_d = {}, {}, {}, {}
        all_dictionaries = {
            "df_main_allGO_hist0": df_class_a,
            "df_main_allAss_hist0": df_class_b,
            "df_main_allCP_GO_hist0": df_class_c,
            "df_main_allCP_Ass_hist0": df_class_d,
        }

        for dict_name, dir_path in directories.items():
            for file_name in os.listdir(dir_path):
                if file_name.endswith(".parquet"):
                    key = file_name.replace(".parquet", "")
                    file_path = os.path.join(dir_path, file_name)
                    try:
                        all_dictionaries[dict_name][key] = pd.read_parquet(
                            file_path, engine="pyarrow"
                        )
                    except Exception as e:
                        print(f"Error reading file {file_path}: {e}")
                        break

        return all_dictionaries

    @staticmethod
    def combine_4_class_dictionaries(
        dict_class_a: Dict[str, pd.DataFrame],
        dict_class_b: Dict[str, pd.DataFrame],
        dict_class_c: Dict[str, pd.DataFrame],
        dict_class_d: Dict[str, pd.DataFrame],
    ) -> pd.DataFrame:
        """
        Combines 4 DataFrame dictionaries into one DataFrame with class labels.
        """

        df_class_a = pd.concat(dict_class_a.values())
        print("Shape of Class A:", df_class_a.shape)
        df_class_b = pd.concat(dict_class_b.values())
        print("Shape of Class B:", df_class_b.shape)
        df_class_c = pd.concat(dict_class_c.values())
        print("Shape of Class C:", df_class_c.shape)
        df_class_d = pd.concat(dict_class_d.values())
        print("Shape of Class D:", df_class_d.shape)

        df_class_a["Class"] = "A"
        df_class_b["Class"] = "B"
        df_class_c["Class"] = "C"
        df_class_d["Class"] = "D"

        df_class_a = df_class_a.drop("-log+b", axis=1)
        df_class_b = df_class_b.drop("-log+b", axis=1)
        df_class_c = df_class_c.drop(
            [
                "GO_TermChild",
                "GO_TermParent",
                "Percentage",
                "Num_Percentage",
                "prelog+b_P",
                "-log+b_P",
            ],
            axis=1,
        )
        df_class_d = df_class_d.drop(
            [
                "GO_TermChild",
                "GO_TermParent",
                "Percentage",
                "Num_Percentage",
                "prelog+b_P",
                "-log+b_P",
            ],
            axis=1,
        )

        combined_df = pd.concat(
            [df_class_a, df_class_b, df_class_c, df_class_d], axis=0, ignore_index=True
        )
        print("Combined DataFrame shape:", combined_df.shape)
        return combined_df

    @staticmethod
    def perform_group_shuffle_split(
        log_handler,
        combined_df: pd.DataFrame,
        output_dir: str,
        n_splits_outer: int = 5,
        n_splits_inner: int = 1,
    ) -> None:
        """
        Performs GroupShuffleSplit on protein names to create train-validation-test splits.

        Results are saved to parquet files under output_dir.
        """
        if log_handler is None:
            log_handler = logging.getLogger(__name__)

        os.makedirs(output_dir, exist_ok=True)

        if "Name" not in combined_df.columns:
            raise ValueError("'Name' column is missing in the DataFrame.")

        groups = combined_df["Name"].values
        outer_splitter = GroupShuffleSplit(
            n_splits=n_splits_outer, test_size=0.2, random_state=42
        )

        for outer_fold, (train_val_idx, test_idx) in enumerate(
            outer_splitter.split(combined_df, groups=groups)
        ):
            train_val_df = combined_df.iloc[train_val_idx].copy()
            test_df = combined_df.iloc[test_idx].copy()

            test_df.to_parquet(
                os.path.join(output_dir, f"CV{outer_fold}_raw_test.parquet")
            )
            train_val_df.to_parquet(
                os.path.join(output_dir, f"CV{outer_fold}_raw_train_val_fold.parquet")
            )

            inner_groups = train_val_df["Name"].values
            inner_splitter = GroupShuffleSplit(
                n_splits=n_splits_inner, test_size=0.2, random_state=42
            )

            for train_idx, val_idx in inner_splitter.split(
                train_val_df, groups=inner_groups
            ):
                train_df = train_val_df.iloc[train_idx].copy()
                val_df = train_val_df.iloc[val_idx].copy()

                log_handler.info("Saved train and val sets for inner split.")

                train_df.to_parquet(
                    os.path.join(output_dir, f"CV{outer_fold}_raw_train.parquet")
                )
                val_df.to_parquet(
                    os.path.join(output_dir, f"CV{outer_fold}_raw_val.parquet")
                )

                log_handler.info("Finished saving inner split sets.")
