#!/bin/bash
#Train: Part5
# Move to the project root (one level up from Scripts_bash)
cd "$(dirname "$0")/.."

# Start timing the entire script
start_time=$(date +%s)

# Define log file locations
log_dir="./Log_files"
mkdir -p "$log_dir"  # Ensure the log directory exists

stderr_file="$log_dir/final_dataset_generator_stderr.txt"

# Print start message
echo "Starting final_dataset_generator training script..."

{
    ###########################################################This part is for Train part6: ONLY#######################################
    ###CV0: Unique Variables per CV
    input1="./Data/Post_CV_Split_Data/train/part5/fold_0_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
    input2="./Data/Post_CV_Split_Data/train/part4/fold_0_train/AA33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/BB33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/CC33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/DD33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/df_AA_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/df_BB_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/df_CC_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/df_DD_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/df_AA_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/df_BB_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/df_CC_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_0_train/df_DD_bin_midpoints.npy"
    input3="./Data/Post_CV_Split_Data/train/part6, fold_0_train,./Log_files,final_dataset_generator_train"
    time python ./Scripts_python/final_dataset_generator.py "$input1" "$input2" "$input3"
    wait
    sleep 30

    ###CV0: Unique Variables per CV
    input1="./Data/Post_CV_Split_Data/train/part5/fold_1_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
    input2="./Data/Post_CV_Split_Data/train/part4/fold_1_train/AA33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/BB33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/CC33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/DD33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/df_AA_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/df_BB_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/df_CC_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/df_DD_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/df_AA_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/df_BB_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/df_CC_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_1_train/df_DD_bin_midpoints.npy"
    input3="./Data/Post_CV_Split_Data/train/part6, fold_1_train,./Log_files,final_dataset_generator_train"
    time python ./Scripts_python/final_dataset_generator.py "$input1" "$input2" "$input3"
    wait
    sleep 30

    ###CV0: Unique Variables per CV
    input1="./Data/Post_CV_Split_Data/train/part5/fold_2_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
    input2="./Data/Post_CV_Split_Data/train/part4/fold_2_train/AA33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/BB33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/CC33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/DD33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/df_AA_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/df_BB_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/df_CC_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/df_DD_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/df_AA_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/df_BB_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/df_CC_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_2_train/df_DD_bin_midpoints.npy"
    input3="./Data/Post_CV_Split_Data/train/part6, fold_2_train,./Log_files,final_dataset_generator_train"
    time python ./Scripts_python/final_dataset_generator.py "$input1" "$input2" "$input3"
    wait
    sleep 30

    ###CV0: Unique Variables per CV
    input1="./Data/Post_CV_Split_Data/train/part5/fold_3_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
    input2="./Data/Post_CV_Split_Data/train/part4/fold_3_train/AA33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/BB33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/CC33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/DD33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/df_AA_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/df_BB_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/df_CC_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/df_DD_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/df_AA_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/df_BB_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/df_CC_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/df_DD_bin_midpoints.npy"
    input3="./Data/Post_CV_Split_Data/train/part6, fold_3_train,./Log_files,final_dataset_generator_train"
    time python ./Scripts_python/final_dataset_generator.py "$input1" "$input2" "$input3"
    wait
    sleep 30

    ###CV0: Unique Variables per CV
    input1="./Data/Post_CV_Split_Data/train/part5/fold_4_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
    input2="./Data/Post_CV_Split_Data/train/part4/fold_4_train/AA33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/BB33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/CC33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/DD33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/df_AA_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/df_BB_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/df_CC_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/df_DD_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/df_AA_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/df_BB_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/df_CC_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train/part4/fold_4_train/df_DD_bin_midpoints.npy"
    input3="./Data/Post_CV_Split_Data/train/part6, fold_4_train,./Log_files,final_dataset_generator_train"
    time python ./Scripts_python/final_dataset_generator.py "$input1" "$input2" "$input3"
    wait
    sleep 30


    ###########################################################This part is for Train_VAL TOGETHER: FOR THE RETRAIN  part6: ONLY#######################################
    ###CV0: Unique Variables per CV
    input1="./Data/Post_CV_Split_Data/train_val/part5/fold_0_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
    input2="./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/AA33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/BB33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/CC33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/DD33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/df_AA_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/df_BB_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/df_CC_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/df_DD_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/df_AA_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/df_BB_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/df_CC_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/df_DD_bin_midpoints.npy"
    input3="./Data/Post_CV_Split_Data/train_val/part6, fold_0_train,./Log_files,final_dataset_generator_train_val"
    time python ./Scripts_python/final_dataset_generator.py "$input1" "$input2" "$input3"
    wait
    sleep 30

    ###CV0: Unique Variables per CV
    input1="./Data/Post_CV_Split_Data/train_val/part5/fold_1_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
    input2="./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/AA33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/BB33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/CC33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/DD33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/df_AA_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/df_BB_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/df_CC_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/df_DD_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/df_AA_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/df_BB_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/df_CC_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/df_DD_bin_midpoints.npy"
    input3="./Data/Post_CV_Split_Data/train_val/part6, fold_1_train,./Log_files,final_dataset_generator_train_val"
    time python ./Scripts_python/final_dataset_generator.py "$input1" "$input2" "$input3"
    wait
    sleep 30

    ###CV0: Unique Variables per CV
    input1="./Data/Post_CV_Split_Data/train_val/part5/fold_2_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
    input2="./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/AA33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/BB33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/CC33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/DD33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/df_AA_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/df_BB_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/df_CC_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/df_DD_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/df_AA_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/df_BB_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/df_CC_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/df_DD_bin_midpoints.npy"
    input3="./Data/Post_CV_Split_Data/train_val/part6, fold_2_train,./Log_files,final_dataset_generator_train_val"
    time python ./Scripts_python/final_dataset_generator.py "$input1" "$input2" "$input3"
    wait
    sleep 30

    ###CV0: Unique Variables per CV
    input1="./Data/Post_CV_Split_Data/train_val/part5/fold_3_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
    input2="./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/AA33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/BB33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/CC33.pkl,\
    ./Data/Post_CV_Split_Data/train/part4/fold_3_train/DD33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/df_AA_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/df_BB_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/df_CC_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/df_DD_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/df_AA_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/df_BB_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/df_CC_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/df_DD_bin_midpoints.npy"
    input3="./Data/Post_CV_Split_Data/train_val/part6, fold_3_train,./Log_files,final_dataset_generator_train_val"
    time python ./Scripts_python/final_dataset_generator.py "$input1" "$input2" "$input3"
    wait
    sleep 30

    ###CV0: Unique Variables per CV
    input1="./Data/Post_CV_Split_Data/train_val/part5/fold_4_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
    input2="./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/AA33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/BB33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/CC33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/DD33.pkl,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/df_AA_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/df_BB_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/df_CC_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/df_DD_bin_edges.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/df_AA_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/df_BB_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/df_CC_bin_midpoints.npy,\
    ./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/df_DD_bin_midpoints.npy"
    input3="./Data/Post_CV_Split_Data/train_val/part6, fold_4_train,./Log_files,final_dataset_generator_train_val"
    time python ./Scripts_python/final_dataset_generator.py "$input1" "$input2" "$input3"
    wait
    sleep 30

}

    # End timing and print execution time to error log
end_time=$(date +%s)
total_time=$((end_time - start_time))
echo "Total Execution Time: $total_time seconds" > "$stderr_file"
echo "Finished final_dataset_generator training script." >> "$stderr_file"












# INPUT1="./JN_gen_final/part5/fold_2_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
# INPUT2="./JN_gen_final/part4/fold_2_train/AA33.pkl,./JN_gen_final/part4/fold_2_train/BB33.pkl,\
# ./JN_gen_final/part4/fold_2_train/CC33.pkl,./JN_gen_final/part4/fold_2_train/DD33.pkl,\
# ./JN_gen_final/part4/fold_2_train/df_AA_bin_edges.npy,./JN_gen_final/part4/fold_2_train/df_BB_bin_edges.npy,\
# ./JN_gen_final/part4/fold_2_train/df_CC_bin_edges.npy,./JN_gen_final/part4/fold_2_train/df_DD_bin_edges.npy,\
# ./JN_gen_final/part4/fold_2_train/df_AA_bin_midpoints.npy,./JN_gen_final/part4/fold_2_train/df_BB_bin_midpoints.npy,\
# ./JN_gen_final/part4/fold_2_train/df_CC_bin_midpoints.npy,./JN_gen_final/part4/fold_2_train/df_DD_bin_midpoints.npy"
# INPUT3="fold_2_train,input_train_data"
# time python python7.py "$INPUT1" "$INPUT2" "$input3"

# INPUT1="./JN_gen_final/part5/fold_3_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
# INPUT2="./JN_gen_final/part4/fold_3_train/AA33.pkl,./JN_gen_final/part4/fold_3_train/BB33.pkl,\
# ./JN_gen_final/part4/fold_3_train/CC33.pkl,./JN_gen_final/part4/fold_3_train/DD33.pkl,\
# ./JN_gen_final/part4/fold_3_train/df_AA_bin_edges.npy,./JN_gen_final/part4/fold_3_train/df_BB_bin_edges.npy,\
# ./JN_gen_final/part4/fold_3_train/df_CC_bin_edges.npy,./JN_gen_final/part4/fold_3_train/df_DD_bin_edges.npy,\
# ./JN_gen_final/part4/fold_3_train/df_AA_bin_midpoints.npy,./JN_gen_final/part4/fold_3_train/df_BB_bin_midpoints.npy,\
# ./JN_gen_final/part4/fold_3_train/df_CC_bin_midpoints.npy,./JN_gen_final/part4/fold_3_train/df_DD_bin_midpoints.npy"
# INPUT3="fold_3_train,input_train_data"
# time python python7.py "$INPUT1" "$INPUT2" "$input3"


# INPUT1="./JN_gen_final/part5/fold_4_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
# INPUT2="./JN_gen_final/part4/fold_4_train/AA33.pkl,./JN_gen_final/part4/fold_4_train/BB33.pkl,\
# ./JN_gen_final/part4/fold_4_train/CC33.pkl,./JN_gen_final/part4/fold_4_train/DD33.pkl,\
# ./JN_gen_final/part4/fold_4_train/df_AA_bin_edges.npy,./JN_gen_final/part4/fold_4_train/df_BB_bin_edges.npy,\
# ./JN_gen_final/part4/fold_4_train/df_CC_bin_edges.npy,./JN_gen_final/part4/fold_4_train/df_DD_bin_edges.npy,\
# ./JN_gen_final/part4/fold_4_train/df_AA_bin_midpoints.npy,./JN_gen_final/part4/fold_4_train/df_BB_bin_midpoints.npy,\
# ./JN_gen_final/part4/fold_4_train/df_CC_bin_midpoints.npy,./JN_gen_final/part4/fold_4_train/df_DD_bin_midpoints.npy"
# INPUT3="fold_4_train,input_train_data"
# time python python7.py "$INPUT1" "$INPUT2" "$input3"

# INPUT1="./JN_gen_final/part5/fold_5_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
# INPUT2="./JN_gen_final/part4/fold_5_train/AA33.pkl,./JN_gen_final/part4/fold_5_train/BB33.pkl,\
# ./JN_gen_final/part4/fold_5_train/CC33.pkl,./JN_gen_final/part4/fold_5_train/DD33.pkl,\
# ./JN_gen_final/part4/fold_5_train/df_AA_bin_edges.npy,./JN_gen_final/part4/fold_5_train/df_BB_bin_edges.npy,\
# ./JN_gen_final/part4/fold_5_train/df_CC_bin_edges.npy,./JN_gen_final/part4/fold_5_train/df_DD_bin_edges.npy,\
# ./JN_gen_final/part4/fold_5_train/df_AA_bin_midpoints.npy,./JN_gen_final/part4/fold_5_train/df_BB_bin_midpoints.npy,\
# ./JN_gen_final/part4/fold_5_train/df_CC_bin_midpoints.npy,./JN_gen_final/part4/fold_5_train/df_DD_bin_midpoints.npy"
# INPUT3="fold_5_train,input_train_data"
# time python python7.py "$INPUT1" "$INPUT2" "$input3"

# INPUT1="./JN_gen_final/part5/fold_hold_out_train/from_part4single_predic_df_B_L_sorted_f_l_with_depth.pkl"
# INPUT2="./JN_gen_final/part4/fold_hold_out_train/AA33.pkl,./JN_gen_final/part4/fold_hold_out_train/BB33.pkl,\
# ./JN_gen_final/part4/fold_hold_out_train/CC33.pkl,./JN_gen_final/part4/fold_hold_out_train/DD33.pkl,\
# ./JN_gen_final/part4/fold_hold_out_train/df_AA_bin_edges.npy,./JN_gen_final/part4/fold_hold_out_train/df_BB_bin_edges.npy,\
# ./JN_gen_final/part4/fold_hold_out_train/df_CC_bin_edges.npy,./JN_gen_final/part4/fold_hold_out_train/df_DD_bin_edges.npy,\
# ./JN_gen_final/part4/fold_hold_out_train/df_AA_bin_midpoints.npy,./JN_gen_final/part4/fold_hold_out_train/df_BB_bin_midpoints.npy,\
# ./JN_gen_final/part4/fold_hold_out_train/df_CC_bin_midpoints.npy,./JN_gen_final/part4/fold_hold_out_train/df_DD_bin_midpoints.npy"
# INPUT3="fold_hold_out_train,input_train_data"
# time python python7.py "$INPUT1" "$INPUT2" "$input3"
