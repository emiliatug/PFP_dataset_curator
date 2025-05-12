#!/bin/bash
#Train: Part5
# Move to the project root (one level up from Scripts_bash)
cd "$(dirname "$0")/.."

# Start timing the entire script
start_time=$(date +%s)

# Define log file locations
log_dir="./Log_files"
mkdir -p "$log_dir"  # Ensure the log directory exists

stderr_file="$log_dir/global1_part456train_stderr.txt"

# Print start message
echo "Starting global1_part456val training script..."
{
    ###########################################################This part is for TEST part: TEST IS RUN ON THE RETRAINED DATA#######################################
    ###CV0: Unique Variables per CV
    ### Define input variables: common variable
    input2="./Data/association.txt,\
    ./Data/category.txt,\
    ./Data/child_parent_percentages,\
    ./Data/GO_definition.txt,\
    ./Data/go_description.txt,\
    ./Data/pygosemsim/_resources/OBO_2022_09_19_go-basic.obo"
    input3="./Data/dataset_protein_seq_df.csv,\
    ./Data/ProteinFeat_df.csv,\
    https://rest.uniprot.org/uniprotkb/"
    input6="./Data/InputGeneratedData/part2p_df_main_allGO_hist0/,\
    ./Data/InputGeneratedData/part2p_df_main_allAss_hist0/,\
    ./Data/InputGeneratedData/part2p_df_main_allCP_GO_hist0/,\
    ./Data/InputGeneratedData/part2p_df_main_allCP_Ass_hist0/"

    ###CV0: Unique Variables per CV
    input1="./Data/InputGeneratedData/input_CV_data/CV0_raw_test.parquet,\
    ./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
    ./Data/pygosemsim/_resources,\
    ./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
    input4="./Data/Post_CV_Split_Data/test/part456test, fold_0_test,./Log_files,global1_part456test, val"
    input5="./Data/Post_CV_Split_Data/train_val/part6/fold_0_train/tmp4_lin_to_analize_df.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part6/fold_0_train/final_train_df.pkl,\
    ../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions"
    time python ./Scripts_python/val_test_processor.py.py "$input1" "$input2" "$input3" "$input4" "$input5", "$input6"
    wait
    sleep 30


    ###CV1: Unique Variables per CV
    input1="./Data/InputGeneratedData/input_CV_data/CV1_raw_test.parquet,\
    ./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
    ./Data/pygosemsim/_resources,\
    ./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
    input4="./Data/Post_CV_Split_Data/test/part456test, fold_1_test,./Log_files,global1_part456test, val"
    input5="./Data/Post_CV_Split_Data/train_val/part6/fold_1_train/tmp4_lin_to_analize_df.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part6/fold_1_train/final_train_df.pkl,\
    ../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions"
    time python ./Scripts_python/val_test_processor.py.py "$input1" "$input2" "$input3" "$input4" "$input5", "$input6"
    wait
    sleep 30

    ###CV2: Unique Variables per CV
    input1="./Data/InputGeneratedData/input_CV_data/CV2_raw_test.parquet,\
    ./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
    ./Data/pygosemsim/_resources,\
    ./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
    input4="./Data/Post_CV_Split_Data/test/part456test, fold_2_test,./Log_files,global1_part456test, val"
    input5="./Data/Post_CV_Split_Data/train_val/part6/fold_2_train/tmp4_lin_to_analize_df.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part6/fold_2_train/final_train_df.pkl,\
    ../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions"
    time python ./Scripts_python/val_test_processor.py.py "$input1" "$input2" "$input3" "$input4" "$input5", "$input6"
    wait
    sleep 30

    ###CV3: Unique Variables per CV
    input1="./Data/InputGeneratedData/input_CV_data/CV3_raw_test.parquet,\
    ./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
    ./Data/pygosemsim/_resources,\
    ./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
    input4="./Data/Post_CV_Split_Data/test/part456test, fold_3_test,./Log_files,global1_part456test, val"
    input5="./Data/Post_CV_Split_Data/train_val/part6/fold_3_train/tmp4_lin_to_analize_df.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part6/fold_3_train/final_train_df.pkl,\
    ../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions"
    time python ./Scripts_python/val_test_processor.py.py "$input1" "$input2" "$input3" "$input4" "$input5", "$input6"
    wait
    sleep 30


    ###CV4: Unique Variables per CV
    input1="./Data/InputGeneratedData/input_CV_data/CV4_raw_test.parquet,\
    ./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
    ./Data/pygosemsim/_resources,\
    ./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
    input4="./Data/Post_CV_Split_Data/test/part456test, fold_4_test,./Log_files,global1_part456test, val"
    input5="./Data/Post_CV_Split_Data/train_val/part6/fold_4_train/tmp4_lin_to_analize_df.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part6/fold_4_train/final_train_df.pkl,\
    ../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions"
    time python ./Scripts_python/val_test_processor.py.py "$input1" "$input2" "$input3" "$input4" "$input5", "$input6"
    wait
    sleep 30

}
# End timing and print execution time to error log
end_time=$(date +%s)
total_time=$((end_time - start_time))
echo "Total Execution Time: $total_time seconds" >> "$stderr_file"
echo "Finished global1_part6 training script." >> "$stderr_file"