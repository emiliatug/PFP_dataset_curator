#!/bin/bash
#Train: Part5
# Move to the project root (one level up from Scripts_bash)
cd "$(dirname "$0")/.."
###########################################################This part is for Train part5: ONLY#######################################
### Define input variables: common variable
input2="./Data/association.txt,\
./Data/category.txt,\
./Data/child_parent_percentages,\
./Data/GO_definition.txt,\
./Data/go_description.txt,\
./Data/pygosemsim/_resources/OBO_2022_09_19_go-basic.obo"
### Define input variables: common variable
input3="./Data/dataset_protein_seq_df.csv,\
./Data/ProteinFeat_df.csv,\
https://rest.uniprot.org/uniprotkb/"

###CV0: Unique Variables per CV
input1="./Data/Post_CV_Split_Data/train/part4/fold_0_train/jonedABCD_PFP.pkl,\
./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
./Data/pygosemsim/_resources,\
./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
input4="./Data/Post_CV_Split_Data/train/part5, fold_0_train,./Log_files,pygosemsim_integrator_train"
time python ./Scripts_python/pygosemsim_integrator.py "$input1" "$input2" "$input3" "$input4"
wait
sleep 30

###CV1: Unique Variables per CV
input1="./Data/Post_CV_Split_Data/train/part4/fold_1_train/jonedABCD_PFP.pkl,\
./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
./Data/pygosemsim/_resources,\
./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
###CV0: Unique Variables per CV
input4="./Data/Post_CV_Split_Data/train/part5, fold_1_train,./Log_files,pygosemsim_integrator_train"
time python ./Scripts_python/pygosemsim_integrator.py "$input1" "$input2" "$input3" "$input4"
wait
sleep 30

###CV2: Unique Variables per CV
input1="./Data/Post_CV_Split_Data/train/part4/fold_2_train/jonedABCD_PFP.pkl,\
./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
./Data/pygosemsim/_resources,\
./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
input4="./Data/Post_CV_Split_Data/train/part5, fold_2_train,./Log_files,pygosemsim_integrator_train"
time python ./Scripts_python/pygosemsim_integrator.py "$input1" "$input2" "$input3" "$input4"
wait
sleep 30

###CV3: Unique Variables per CV
input1="./Data/Post_CV_Split_Data/train/part4/fold_3_train/jonedABCD_PFP.pkl,\
./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
./Data/pygosemsim/_resources,\
./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
input4="./Data/Post_CV_Split_Data/train/part5, fold_3_train,./Log_files,pygosemsim_integrator_train"
time python ./Scripts_python/pygosemsim_integrator.py "$input1" "$input2" "$input3" "$input4"
wait
sleep 30

###CV4: Unique Variables per CV
input1="./Data/Post_CV_Split_Data/train/part4/fold_4_train/jonedABCD_PFP.pkl,\
./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
./Data/pygosemsim/_resources,\
./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
input4="./Data/Post_CV_Split_Data/train/part5, fold_4_train,./Log_files,pygosemsim_integrator_train"
time python ./Scripts_python/pygosemsim_integrator.py "$input1" "$input2" "$input3" "$input4"
wait
sleep 30

###########################################################This part is for Train_VAL TOGETHER: FOR THE RETRAIN  part5: ONLY#######################################
### Define input variables: common variable
input2="./Data/association.txt,\
./Data/category.txt,\
./Data/child_parent_percentages,\
./Data/GO_definition.txt,\
./Data/go_description.txt,\
./Data/pygosemsim/_resources/OBO_2022_09_19_go-basic.obo"
### Define input variables: common variable
input3="./Data/dataset_protein_seq_df.csv,\
./Data/ProteinFeat_df.csv,\
https://rest.uniprot.org/uniprotkb/"

###CV0: Unique Variables per CV
input1="./Data/Post_CV_Split_Data/train_val/part4/fold_0_train/jonedABCD_PFP.pkl,\
./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
./Data/pygosemsim/_resources,\
./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
input4="./Data/Post_CV_Split_Data/train_val/part5, fold_0_train,./Log_files,pygosemsim_integrator_train_val"
time python ./Scripts_python/pygosemsim_integrator.py "$input1" "$input2" "$input3" "$input4"
wait
sleep 30

###CV1: Unique Variables per CV
input1="./Data/Post_CV_Split_Data/train_val/part4/fold_1_train/jonedABCD_PFP.pkl,\
./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
./Data/pygosemsim/_resources,\
./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
###CV0: Unique Variables per CV
input4="./Data/Post_CV_Split_Data/train_val/part5, fold_1_train,./Log_files,pygosemsim_integrator_train_val"
time python ./Scripts_python/pygosemsim_integrator.py "$input1" "$input2" "$input3" "$input4"
wait
sleep 30

###CV2: Unique Variables per CV
input1="./Data/Post_CV_Split_Data/train_val/part4/fold_2_train/jonedABCD_PFP.pkl,\
./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
./Data/pygosemsim/_resources,\
./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
input4="./Data/Post_CV_Split_Data/train_val/part5, fold_2_train,./Log_files,pygosemsim_integrator_train_val"
time python ./Scripts_python/pygosemsim_integrator.py "$input1" "$input2" "$input3" "$input4"
wait
sleep 30

###CV3: Unique Variables per CV
input1="./Data/Post_CV_Split_Data/train_val/part4/fold_3_train/jonedABCD_PFP.pkl,\
./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
./Data/pygosemsim/_resources,\
./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
input4="./Data/Post_CV_Split_Data/train_val/part5, fold_3_train,./Log_files,pygosemsim_integrator_train_val"
time python ./Scripts_python/pygosemsim_integrator.py "$input1" "$input2" "$input3" "$input4"
wait
sleep 30

###CV4: Unique Variables per CV
input1="./Data/Post_CV_Split_Data/train_val/part4/fold_4_train/jonedABCD_PFP.pkl,\
./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,\
./Data/pygosemsim/_resources,\
./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211"
input4="./Data/Post_CV_Split_Data/train_val/part5, fold_4_train,./Log_files,pygosemsim_integrator_train_val"
time python ./Scripts_python/pygosemsim_integrator.py "$input1" "$input2" "$input3" "$input4"
wait
sleep 30





#!/bin/bash

#first 5 calls are 5 CV and the last call is the full one

# time python python5.py './Data/Post_CV_Split_Data/train/part4/jonedABCD_PFP.pkl,./Data/updated_uniprot_6_species_with_GOterms_GAF.csv,./Data/pygosemsim/_resources,./Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211' './Data/association.txt,./Data/category.txt,./Data/child_parent_percentages,./Data/GO_definition.txt,./Data/go_description.txt,./Data/pygosemsim/_resources/OBO_2022_09_19_go-basic.obo' './Data/dataset_protein_seq_df.csv,./Data/ProteinFeat_df.csv,https://rest.uniprot.org/uniprotkb/' 'fold_0_train','./Log_files,global1_part4' 

# time python python5.py './JN_gen_final/part4/fold_2_train/jonedABCD_PFP.pkl,../Data/updated_uniprot_6_species_with_GOterms_GAF.csv,../Data/pygosemsim/_resources,../Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211' '../Data/association.txt,../Data/category.txt,../Data/child_parent_percentages,../Data/GO_definition.txt,../Data/go_description.txt,./Datapygosemsim/_resources/OBO_2022_09_19_go-basic.obo' '../Data/dataset_protein_seq_df.csv,../Data/ProteinFeat_df.csv,https://rest.uniprot.org/uniprotkb/' 'fold_2_train'

# time python python5.py './JN_gen_final/part4/fold_3_train/jonedABCD_PFP.pkl,../Data/updated_uniprot_6_species_with_GOterms_GAF.csv,../Data/pygosemsim/_resources,../Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211' '../Data/association.txt,../Data/category.txt,../Data/child_parent_percentages,../Data/GO_definition.txt,../Data/go_description.txt,./Datapygosemsim/_resources/OBO_2022_09_19_go-basic.obo' '../Data/dataset_protein_seq_df.csv,../Data/ProteinFeat_df.csv,https://rest.uniprot.org/uniprotkb/' 'fold_3_train'

# time python python5.py './JN_gen_final/part4/fold_4_train/jonedABCD_PFP.pkl,../Data/updated_uniprot_6_species_with_GOterms_GAF.csv,../Data/pygosemsim/_resources,../Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211' '../Data/association.txt,../Data/category.txt,../Data/child_parent_percentages,../Data/GO_definition.txt,../Data/go_description.txt,./Datapygosemsim/_resources/OBO_2022_09_19_go-basic.obo' '../Data/dataset_protein_seq_df.csv,../Data/ProteinFeat_df.csv,https://rest.uniprot.org/uniprotkb/' 'fold_4_train'

# time python python5.py './JN_gen_final/part4/fold_5_train/jonedABCD_PFP.pkl,../Data/updated_uniprot_6_species_with_GOterms_GAF.csv,../Data/pygosemsim/_resources,../Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211' '../Data/association.txt,../Data/category.txt,../Data/child_parent_percentages,../Data/GO_definition.txt,../Data/go_description.txt,./Datapygosemsim/_resources/OBO_2022_09_19_go-basic.obo' '../Data/dataset_protein_seq_df.csv,../Data/ProteinFeat_df.csv,https://rest.uniprot.org/uniprotkb/' 'fold_5_train'

# time python python5.py './JN_gen_final/part4/fold_hold_out_train/jonedABCD_PFP.pkl,../Data/updated_uniprot_6_species_with_GOterms_GAF.csv,../Data/pygosemsim/_resources,../Data/pygosemsim/_resources/cleaned_noTrEMBL_two_col_goa_uniprot_all.gpi.211' '../Data/association.txt,../Data/category.txt,../Data/child_parent_percentages,../Data/GO_definition.txt,../Data/go_description.txt,./Datapygosemsim/_resources/OBO_2022_09_19_go-basic.obo' '../Data/dataset_protein_seq_df.csv,../Data/ProteinFeat_df.csv,https://rest.uniprot.org/uniprotkb/' 'fold_hold_out_train'


# ./Scripts_bash/run_global1_part4.sh