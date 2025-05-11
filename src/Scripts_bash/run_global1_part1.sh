#!/bin/bash

# Move to the project root (one level up from Scripts_bash)
cd "$(dirname "$0")/.."

# Define input variables
input1="../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions,\
../../../07122024DSets_Tr/DataSet_Evalues_CSV/ProteinEvaluesGOterms,\
../../../07122024_Tr/PFP*/tmp/"

input2="https://rest.uniprot.org/uniprotkb/"

input3="./Data/InputGeneratedData/part1p_df_ProtEvalGO_upd/,\
./Data/InputGeneratedData/part1p_df_ProtEvalGO_updGO/,\
./Data/InputGeneratedData/part1p_df_ProtEvalGO_updAss/"

input4="../../../07122024DSets_Tr/DataSet_Evalues_CSV/Parent_dir,\
../../../07122024_Tr/PFP1/bin/data_pfp/child_parent_percentages,\
./Data/InputGeneratedData/part1p_df_child_parent/"

input5="./Data/InputGeneratedData/part1p_df_main_allGO/,\
./Data/InputGeneratedData/part1p_df_main_allAss/,\
./Data/InputGeneratedData/part1p_df_main_allCP_GO/,\
./Data/InputGeneratedData/part1p_df_main_allCP_Ass/"

input6="./Log_files,global1_part1"

# Run the Python script
time python ./Scripts_python/class_generator.py "$input1" "$input2" "$input3" "$input4" "$input5" "$input6"