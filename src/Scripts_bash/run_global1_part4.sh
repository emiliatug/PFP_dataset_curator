#!/bin/bash
cd "$(dirname "$0")/.."
###########################################################This part is for Train part4: ONLY#######################################
### Define input variables: common variable
input1="./Data/InputGeneratedData/part2p_df_main_allGO_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allAss_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allCP_GO_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allCP_Ass_hist0/"

###CV0: Unique Variables per CV
input2="./Data/InputGeneratedData/input_CV_data/CV0_raw_train.parquet,../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions" 
input3="./Data/Post_CV_Split_Data/train/part4,fold_0_train,./Log_files,create_sparse_training_arrays_train"
time python ./Scripts_python/create_sparse_training_arrays.py "$input1" "$input2" "$input3"
wait
sleep 30

##CV1: Unique Variables per CV
input2="./Data/InputGeneratedData/input_CV_data/CV1_raw_train.parquet,../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions" 
input3="./Data/Post_CV_Split_Data/train/part4,fold_1_train,./Log_files,create_sparse_training_arrays_train"
time python ./Scripts_python/create_sparse_training_arrays.py "$input1" "$input2" "$input3"
wait
sleep 30

###CV2: Unique Variables per CV
input2="./Data/InputGeneratedData/input_CV_data/CV2_raw_train.parquet,../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions" 
input3="./Data/Post_CV_Split_Data/train/part4,fold_2_train,./Log_files,create_sparse_training_arrays_train"
time python ./Scripts_python/create_sparse_training_arrays.py "$input1" "$input2" "$input3"
wait
sleep 30

###CV3: Unique Variables per CV
input2="./Data/InputGeneratedData/input_CV_data/CV3_raw_train.parquet,../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions" 
input3="./Data/Post_CV_Split_Data/train/part4,fold_3_train,./Log_files,create_sparse_training_arrays_train"
time python ./Scripts_python/create_sparse_training_arrays.py "$input1" "$input2" "$input3"
wait
sleep 30

###CV4: Unique Variables per CV
input2="./Data/InputGeneratedData/input_CV_data/CV4_raw_train.parquet,../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions" 
input3="./Data/Post_CV_Split_Data/train/part4,fold_4_train,./Log_files,create_sparse_training_arrays_train"
time python ./Scripts_python/create_sparse_training_arrays.py "$input1" "$input2" "$input3"
wait
sleep 30


###########################################################This part is for Train_VAL TOGETHER: FOR THE RETRAIN  part4: ONLY#######################################
### Define input variables: common variable
input1="./Data/InputGeneratedData/part2p_df_main_allGO_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allAss_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allCP_GO_hist0/,\
./Data/InputGeneratedData/part2p_df_main_allCP_Ass_hist0/"

###CV0: Unique Variables per CV
input2="./Data/InputGeneratedData/input_CV_data/CV0_raw_train_val_fold.parquet,../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions" 
input3="./Data/Post_CV_Split_Data/train_val/part4,fold_0_train,./Log_files,create_sparse_training_arrays_train_val"
time python ./Scripts_python/create_sparse_training_arrays.py "$input1" "$input2" "$input3"
wait
sleep 30

###CV1: Unique Variables per CV
input2="./Data/InputGeneratedData/input_CV_data/CV1_raw_train_val_fold.parquet,../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions" 
input3="./Data/Post_CV_Split_Data/train_val/part4,fold_1_train,./Log_files,create_sparse_training_arrays_train_val"
time python ./Scripts_python/create_sparse_training_arrays.py "$input1" "$input2" "$input3"
wait
sleep 30

###CV2: Unique Variables per CV
input2="./Data/InputGeneratedData/input_CV_data/CV2_raw_train_val_fold.parquet,../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions" 
input3="./Data/Post_CV_Split_Data/train_val/part4,fold_2_train,./Log_files,create_sparse_training_arrays_train_val"
time python ./Scripts_python/create_sparse_training_arrays.py "$input1" "$input2" "$input3"
wait
sleep 30

###CV3: Unique Variables per CV
input2="./Data/InputGeneratedData/input_CV_data/CV3_raw_train_val_fold.parquet,../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions" 
input3="./Data/Post_CV_Split_Data/train_val/part4,fold_3_train,./Log_files,create_sparse_training_arrays_train_val"
time python ./Scripts_python/create_sparse_training_arrays.py "$input1" "$input2" "$input3"
wait
sleep 30

###CV4: Unique Variables per CV
input2="./Data/InputGeneratedData/input_CV_data/CV4_raw_train_val_fold.parquet,../../../07122024DSets_Tr/DataSet_Evalues_CSV/GO_predictions" 
input3="./Data/Post_CV_Split_Data/train_val/part4,fold_4_train,./Log_files,create_sparse_training_arrays_train_val"
time python ./Scripts_python/create_sparse_training_arrays.py "$input1" "$input2" "$input3"
wait
sleep 30
