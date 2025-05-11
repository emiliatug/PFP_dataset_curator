import itertools
from typing import Any, Dict, Tuple, Hashable

import numpy as np
from pandas import DataFrame
from scipy.sparse import csr_matrix, vstack
from tqdm import tqdm

from Utils.utils import *


class SparseTrainingArrayHandler:
    """
    Class to handle sparse training array generation and related operations.
    """

    @staticmethod
    def split_dataframe_by_class(
        df_input: pd.DataFrame,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Splits a DataFrame into four DataFrames based on 'Class' column values A, B, C, D.
        """
        df_a = df_input[df_input["Class"] == "A"]
        df_b = df_input[df_input["Class"] == "B"]
        df_c = df_input[df_input["Class"] == "C"]
        df_d = df_input[df_input["Class"] == "D"]
        return df_a, df_b, df_c, df_d

    @staticmethod
    def enrich_dataframe_by_class(
        df_input: pd.DataFrame, full_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Checks for leakage by ensuring only matching 'Name' values are kept from full_df.
        """
        names = df_input["Name"].drop_duplicates().to_list()
        return full_df[full_df["Name"].isin(names)]

    class RangeCalculator:
        def __init__(self) -> None:
            self.call_count = 1

        def give_range(
            self, df: pd.DataFrame, col_name: str
        ) -> Tuple[List[Any], List[Any]]:
            ranges = [[df[col_name].min(), df[col_name].max()]]
            evalue_range = [
                [
                    df[df[col_name] == df[col_name].min()]["Num_found_evalue"]
                    .drop_duplicates()
                    .iloc[0],
                    df[df[col_name] == df[col_name].max()]["Num_found_evalue"]
                    .drop_duplicates()
                    .iloc[0],
                ]
            ]
            self.call_count += 1
            return ranges[0], evalue_range[0]

        def reset_count(self) -> None:
            self.call_count = 1

        @staticmethod
        def find_min_max(list1, list2, list3, list4) -> List[Any]:
            values = list(itertools.chain(list1, list2, list3, list4))
            return [min(values), max(values)]

    @staticmethod
    def make_train_df_to_dic(
        df_input: pd.DataFrame,
    ) -> dict[Hashable, DataFrame]:
        """
        Groups dataframe by 'Name' and returns dictionary of sub-dataframes.
        """

        return_dic = {
            name: group.reset_index(drop=True)
            for name, group in df_input.groupby("Name")
        }

        return return_dic

    @staticmethod
    def get_histogram_values_sparse(
        input_dict: Dict[str, pd.DataFrame],
        value_col: str,
        binning_method: str,
        range_limits: List[float],
        suffix: str,
    ) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
        """
        Generates histogram-based sparse array for each protein grouped by GO term.
        """
        all_values = np.concatenate(
            [df[value_col].values.flatten() for df in input_dict.values()]
        )
        bin_edges = np.histogram_bin_edges(
            all_values, bins=binning_method, range=range_limits
        )
        bin_midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2

        result = {}
        for count, (key, df) in enumerate(input_dict.items()):
            grouped = df.groupby("GO_Term")[value_col].apply(list)
            data = {"Name": [], "GO": [], f"filtered_count{suffix}": []}
            for go_term, values in grouped.items():
                counts = np.histogram(values, bins=bin_edges)[0]
                sparse_counts = pd.arrays.SparseArray(counts)
                data["Name"].append(key)
                data["GO"].append(go_term)
                data[f"filtered_count{suffix}"].append(sparse_counts)
            result[count] = pd.DataFrame(data)

        final_df = pd.concat(result.values())
        return final_df, bin_edges, bin_midpoints

    @staticmethod
    def sparse_identifier(array: Any) -> Tuple[Any, ...]:
        """
        Creates a hashable identifier for a SparseArray.
        """
        if isinstance(array, pd.arrays.SparseArray):
            return tuple((i, array[i]) for i in array.sp_index.to_int_index().indices)
        return tuple(array)

    def drop_duplicates_sparse(
        self, df_input: pd.DataFrame, sparse_col: str
    ) -> pd.DataFrame:
        """
        Drops duplicates based on SparseArray content.
        """
        df_input["TT"] = df_input[sparse_col].apply(self.sparse_identifier)
        df_input = df_input.drop_duplicates(subset=["Name", "GO", "TT"])
        return df_input.drop(columns=["TT"])

    @staticmethod
    def convert_sparse_column_to_dataframe(
        df_input: pd.DataFrame, sparse_col: str
    ) -> pd.DataFrame:
        """
        Converts a column of SparseArrays to a dense sparse matrix DataFrame.
        """
        sparse_matrices = [
            csr_matrix(arr)
            for arr in tqdm(df_input[sparse_col], desc="Processing Sparse Arrays")
        ]
        combined_matrix = vstack(sparse_matrices)
        combined_csr = csr_matrix(combined_matrix)
        sparse_df = pd.DataFrame.sparse.from_spmatrix(combined_csr)
        return pd.concat(
            [df_input.reset_index(drop=True).drop(columns=[sparse_col]), sparse_df],
            axis=1,
        )

    @staticmethod
    def merge_sparse_columns_fixed(
        df_input: pd.DataFrame, target_cols: int
    ) -> pd.DataFrame:
        """
        Merges lowest pairs of columns until target column count is reached.
        """
        numeric_data = df_input._get_numeric_data()
        numeric_data = numeric_data.apply(
            lambda x: pd.arrays.SparseArray(x, fill_value=0)
        )

        while numeric_data.shape[1] > target_cols:
            sums = [numeric_data[col].sum() for col in numeric_data.columns]
            pair_sums = np.array(sums[:-1]) + np.array(sums[1:])
            min_idx = int(np.argmin(pair_sums))
            col1, col2 = (
                numeric_data.columns[min_idx],
                numeric_data.columns[min_idx + 1],
            )
            merged = pd.arrays.SparseArray(
                numeric_data[col1] + numeric_data[col2], fill_value=0
            )
            merged_name = f"M_{col1}_{col2}"
            numeric_data.insert(min_idx, merged_name, merged)
            numeric_data.drop([col1, col2], axis=1, inplace=True)

        return pd.concat([df_input.iloc[:, :2], numeric_data], axis=1)

    @staticmethod
    def add_sum_row(df_input: pd.DataFrame) -> pd.DataFrame:
        """
        Adds a row of summed values for numeric columns.
        """
        numeric_cols = df_input.select_dtypes(include=[np.number]).columns
        sums = pd.arrays.SparseArray(
            df_input[numeric_cols].sum(axis=0).astype("int64"), fill_value=0
        )
        sums_df = pd.DataFrame([sums], columns=numeric_cols, index=["sum"])
        return pd.concat([df_input, sums_df])

    @staticmethod
    def merge_sparse_columns_temp(
        df_input: pd.DataFrame, target_cols: int
    ) -> pd.DataFrame:
        """
        Alternative merge function using SparseArrays directly.
        """
        numeric_data = df_input.apply(lambda x: pd.arrays.SparseArray(x, fill_value=0))

        while numeric_data.shape[1] > target_cols:
            sums = [numeric_data[col].sum() for col in numeric_data.columns]
            pair_sums = np.array(sums[:-1]) + np.array(sums[1:])
            min_idx = int(np.argmin(pair_sums))
            col1, col2 = (
                numeric_data.columns[min_idx],
                numeric_data.columns[min_idx + 1],
            )
            merged = pd.arrays.SparseArray(
                numeric_data[col1] + numeric_data[col2], fill_value=0
            )
            merged_name = f"M_{col1}_{col2}"
            numeric_data.insert(min_idx, merged_name, merged)
            numeric_data.drop([col1, col2], axis=1, inplace=True)

        return numeric_data

    def group_sparse_matrix(
        self, df_input: pd.DataFrame, batch_size: int
    ) -> pd.DataFrame:
        """
        Groups a sparse matrix dataframe into smaller batches and merges pairs until final structure is stable.
        """
        columns = list(df_input.columns)
        result = {}
        for i in range(0, len(columns), batch_size):
            batch_df = df_input[columns[i : i + batch_size]]
            merged = self.merge_sparse_columns_temp(batch_df, target_cols=100)
            merged_with_sum = self.add_sum_row(merged)
            result[i // batch_size] = merged_with_sum

        final_df = pd.concat(result.values(), axis=1)

        if final_df.shape[1] == 100:
            return final_df
        else:
            return self.group_sparse_matrix(final_df, batch_size=400)

    @staticmethod
    def generate_original_predictions(directory: str) -> Dict[str, pd.DataFrame]:
        """
        Reads PFP original prediction files and filters by Raw_Score > 100.
        """
        result = {}
        for filename in os.listdir(directory):
            if filename.endswith(".csv"):
                filepath = os.path.join(directory, filename)
                df = pd.read_csv(
                    filepath,
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
                filtered = df[df["Raw_Score"] > 100.00]
                if not filtered.empty:
                    name = filename[:-4]
                    filtered["Name"] = name
                    result[name] = filtered.copy()

        return result
