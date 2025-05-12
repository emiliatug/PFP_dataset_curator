import pandas as pd

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)

from Utils.utils import *


class ValTestCountHandler:
    def __init__(self):
        pass

    @staticmethod
    def adding_raw_counts_to_bins(input_df, reference_train_input, var1, var2):
        """This fx takes the test or val input and organizes the row (by gettings its counts), according to the train's data bins"""
        filtered_df = input_df[
            (input_df["GO_Term"] == "0055093") & (input_df["Name"] == "Q86C56")
        ]

        col_l = []
        column_names_list = (
            reference_train_input.columns.to_list()
        )  # length of 357, a reference for all the classes
        column_names_list1 = [
            x for x in column_names_list if var1 in x
        ]  # just select bins for that specific class
        column_names_list2 = [
            x[:-1] for x in column_names_list1
        ]  # cleans the name of the bin, by removing the class association
        column_names_list2 = [(float(a), float(b)) for a, b in column_names_list2]
        bin_edges, new_columns = [], []

        for min_score in column_names_list2:
            for el in min_score:
                bin_edges.append(float(el))
        bin_edges[-1] += 1e-10
        bin_edges = sorted(bin_edges)
        input_df_tmp = input_df.copy()
        input_df_tmp["Ranges"] = pd.cut(
            input_df[var2],
            bins=bin_edges,
            labels=column_names_list2,
            right=False,
            duplicates="drop",
        )
        input_df_tmp = input_df_tmp.dropna()

        # for the big_df: include the values that are bigger than the right bound in the last bin AND equal to the right bound
        big_tmp = input_df[input_df[var2].to_numpy() >= bin_edges[-1]]
        big_tmp = big_tmp.copy()
        big_tmp["Ranges"] = [tuple(bin_edges[-2:])] * len(big_tmp)

        # for the small_df: only include values that are smaller than the elft bound for the smallest bin
        small_tmp = input_df[input_df[var2].to_numpy() < bin_edges[0]]
        small_tmp = small_tmp.copy()
        small_tmp["Ranges"] = [tuple(bin_edges[:2])] * len(small_tmp)
        tmp1 = pd.concat(
            [df for df in [input_df_tmp, small_tmp, big_tmp] if not df.empty], axis=0
        ).reset_index(drop=True)
        tmp1 = (
            tmp1.groupby(
                ["GO_Term", "Ontology", "Raw_Score", "Rank", "Name", "Ranges"],
                observed=True,
            )
            .size()
            .reset_index(name="Count")
        )
        tmp1 = tmp1.pivot(
            index=["GO_Term", "Ontology", "Raw_Score", "Rank", "Name"],
            columns="Ranges",
            values="Count",
        ).fillna(0)
        tmp1 = tmp1.reindex(columns=column_names_list2, fill_value=0).reset_index()

        for idx, col in enumerate(tmp1.columns):
            if isinstance(col, tuple):
                transformed_col = (
                    "Col:[" + str(col[0]) + "-" + str(col[1]) + "]:" + var1
                )
                new_columns.append(transformed_col)
            else:
                new_columns.append(col)
        tmp1.columns = new_columns
        tmp1 = tmp1.reset_index(drop=True)
        return tmp1

    @staticmethod
    def filter_out_root_terms(input_df, col_list, list_root):
        """This fx filters out the root terms for multiple columns"""
        for col_name in col_list:
            input_df = input_df[~input_df[col_name].isin(list_root)]
        return input_df

    # filtered_df=filter_out_root_terms(input_df, col_list,list_root)

    def clean_val_or_df(self, input1, refer_df, col_nam1, list_root):
        """This fx cleans the validation or test df to be like train df by a)changing the column name b) removing extra columns c)reordering
        columns. The reference df is the final training df: Using it only for the Col's name, and order of Columns. There is no Leakage!
        """

        Count = 0
        col_name_list, col_to_drop = [], []

        # part1: Renames columns on input_df
        for col_name in input1.columns:
            if "Col:[" in col_name and "]:" in col_name:
                col_name_list.append(f"Feat_Freq{Count}")
                Count = Count + 1
            else:
                col_name_list.append(col_name)
        input1.columns = col_name_list
        # part2: Drops columns that are not in the training df
        for col_name in input1.columns:
            if col_name not in refer_df.columns:
                col_to_drop.append(col_name)
        input1 = input1.drop(col_to_drop, axis=1)
        # part3: CCalling "filter_out_root_terms" to get rid of any roots if they are presetn
        filtered_df = self.filter_out_root_terms(input1, col_nam1, list_root)

        return_df = input1[refer_df.columns]

        return return_df
