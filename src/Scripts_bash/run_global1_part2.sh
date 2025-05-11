#!/bin/bash

# Move to the project root (one level up from Scripts_bash)
cd "$(dirname "$0")/.."

# Define input variables
input1="./Data/InputGeneratedData/part1p_df_main_allGO/,\
./Data/InputGeneratedData/part1p_df_main_allAss/,\
./Data/InputGeneratedData/part1p_df_main_allCP_GO/,\
./Data/InputGeneratedData/part1p_df_main_allCP_Ass/,\
./Data/InputGeneratedData/part1p_df_child_parent/"

input2="./Data/InputGeneratedData/part2p_df_main_allGO_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allAss_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allCP_GO_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allCP_Ass_hist0/"

input3="./Log_files",\
"evalue_generator"
# Run the Python script with the inputs
time python ./Scripts_python/evalue_generator.py "$input1" "$input2" "$input3"