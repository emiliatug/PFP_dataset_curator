import itertools
import re

import numpy as np
import pandas as pd

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)

import os
import sys

# Now import modules
from Utils.utils import *


class FinalDatasetGenerationHandler:
    def __init__(self):
        pass

    ###Fx from the JN7
    @staticmethod
    def combine_scores_and_counts(input_df, df2):
        """This takes the counts (of 100 features of raw daa) and combines with the previously generated .depth dataframe"""

        input_df = input_df.rename(
            columns={"Protein": "Name", "Original_PFP_Pred_GO_term": "GO"}
        )
        return_df = pd.merge(input_df, df2, on=["Name", "GO"], how="inner")
        return return_df

    @staticmethod
    def find_midpoints(l1):
        """This just takes the column names and organizes them in the list of lists: [0,1,[2,3]...], where 1st two columns are 0,1 and third
        is merged bins 2,3"""

        midpoint_l = []
        for unit in l1[2:]:
            if isinstance(unit, int):
                midpoint_l.append(unit)
            if isinstance(unit, str):
                if "M" in unit:
                    numbers = re.findall(r"\d+", unit)
                    midpoint_l.append(list(map(int, numbers)))
        return midpoint_l

    @staticmethod
    def read_clean_part4_var(input_df):
        input_df = input_df.dropna(subset=["Name", "GO"]).copy()
        input_df["GO"] = "GO:" + input_df["GO"]
        return input_df

    @staticmethod
    def make_ranges_from_bins(index_list, bin_edges):
        """There are 3 returns for this function:new_column_ranges is the list of all column names(list of lists) except the bins numbers
        are replaced with the actual values (-logEvalues), while,return_list are the borders of each bin. Given a list of lists,
        the left most border is left, while the right most border is right. Now,the return_list,updated_2 is used for graphing visuals.
        """
        new_column_ranges = []
        for x in index_list:
            if isinstance(x, int):
                new_column_ranges.append([bin_edges[x], bin_edges[x + 1]])

            if isinstance(x, list):
                for i, group in enumerate(index_list):
                    if isinstance(group, int):
                        index_list[i] = [group]
                tmpL = []
                for y in x:
                    if y == index_list[-1][-1]:
                        tmpL.append([bin_edges[y], bin_edges[y + 1]])
                        break
                    tmpL.append([bin_edges[y], bin_edges[y + 1]])
                new_column_ranges.append(tmpL)

        return_list = []
        for x in new_column_ranges:
            if isinstance(x[0], float):
                x = list(np.around(np.array(x), 5))
                return_list.append(x)

            if isinstance(x[0], list):
                tmp_list = list(itertools.chain.from_iterable(x))

                tmp_list1 = [tmp_list[0], tmp_list[-1]]
                tmp_list2 = list(np.around(np.array(tmp_list1), 5))
                return_list.append(tmp_list2)

        return new_column_ranges, return_list

    ##updated_2 is the list that I generated for the plotting (start, extend by), (start, extend by)
    @staticmethod
    def rename_df(input1, input2):
        together1 = []
        return_df = input1.copy(deep=True)
        first = list(input1.columns[:7])
        together = first + input2
        for x in together:
            if type(x) == list:
                together1.append(tuple(x))
            else:
                together1.append(x)
        return_df.columns = together1

        return_df = return_df.drop(["BinaryScore"], axis=1)
        return return_df

    @staticmethod
    def make_ranges_for_bgraph(input_list):
        """This fx is just making bins to be able to graph them below"""

        updated_l = []
        updated_2 = []
        for x in input_list:
            updated_l.append([x[0], x[1] - x[0]])
            updated_l = list(np.around(np.array(updated_l), 5))
        for x in updated_l:
            updated_2.append(tuple(x))

        return updated_2

    @staticmethod
    def plot_ranges(ax, ranges, color_map, label):
        for i, (start, length) in enumerate(ranges):
            ax.broken_barh(
                [(start, length)],
                (0.1, 0.8),
                facecolors=color_map(i),
                edgecolor="black",
                linewidth=1,
            )
        ax.set_xlim(9e-05, 182.097)  # Adjust based on your data
        ax.set_xlabel("Number Line")
        ax.set_title(label)
        ax.get_yaxis().set_visible(False)
