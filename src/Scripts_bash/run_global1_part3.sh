#!/bin/bash

# Move to the project root (one level up from Scripts_bash)
cd "$(dirname "$0")/.."

# Define input variables
input1="./Data/InputGeneratedData/part2p_df_main_allGO_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allAss_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allCP_GO_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allCP_Ass_hist0/"

input2="./Data/InputGeneratedData/input_CV_data"
input3="./Log_files,global1_part3"


# Run the Python script
time python ./Scripts_python/combined_train_val_generator.py "$input1" "$input2" "$input3"