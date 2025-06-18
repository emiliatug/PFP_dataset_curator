# Protein Function Prediction (PFP) Dataset Curation 

## Background
This project (Protein Function Prediction Dataset Curation and pfp_prediction_evaluator) implements the full pipeline for generating, curating, and labeling protein dataset to evaluate the performance of the Protein Function Prediction (PFP) method developed by the Kihara Lab, (https://kiharalab.org/). 

While traditionally protein functions have been determined experimentally. Since 1990s there has been a lot of progress in computational protein function prediction research. Normally the predictions are either evaluated on benchmark dataset, or through challenges, like [CAFA](https://biofunctionprediction.org/cafa/), where the organizers provide sequence data and different teams submit their predictions which are validated by organizers through experimental annotations.

This page (Protein Function Prediction (PFP) Dataset Curation) contains a series of pipelines for data preprocessing and curating for 9001 proteins which I am using in the [pfp_prediction_evaluator](https://github.com/emiliatug/pfp_prediction_evaluator), to assess the accuracy of the PFP method through machine learning. 

### Data Generation and Curation
This repository focuses on building such a curated dataset through:
protein and Gene Ontology Terms from [UniProt](https://www.uniprot.org/), and [Gene Ontology Consortium](https://geneontology.org/) feature extraction for the list of pre-selected 9001 proteins
Labeling each prediction using the Lin Semantic Similarity Score

The cluster-representative proteins were selected using sequences derived from the [UniClust](https://uniclust.mmseqs.com/) database, based on 30% sequence identity and filtered to only include SwissProt proteins. These proteins serve as a representative subset of UniProt’s protein space and were generated in bash. 

### Structured Dataset
For each protein–PFP-predicted GO term pair, I assembled:

-The frequency of the PSI-BLAST hits from the PFP's intermediate outputs

-Supporting metadata from both UniProt and GO consortium in either one hot-encodings or numerical features

-A Lin score label comparing the predicted GO term to the correct set of annotations

## Environment

Please activate the environment before starting:

```bash
pip install requirements.txt
```

---

## Execution Order

### Step 1: run\_global1\_part1.sh

```bash
chmod +x Scripts_bash/run_global1_part1.sh
./Scripts_bash/run_global1_part1.sh
```

### Step 2: run\_global1\_part2.sh 

```bash
chmod +x Scripts_bash/run_global1_part2.sh
./Scripts_bash/run_global1_part2.sh
```

### Step 3: run\_global1\_part3.sh 

```bash
chmod +x Scripts_bash/run_global1_part3.sh
./Scripts_bash/run_global1_part3.sh
```

### Step 4: run\_global1\_part4.sh 

```bash
chmod +x Scripts_bash/run_global1_part4.sh
./Scripts_bash/run_global1_part4.sh
```

### Step 5: run\_global1\_part5.sh 

```bash
chmod +x Scripts_bash/run_global1_part5.sh
./Scripts_bash/run_global1_part5.sh
```

### Step 6: run\_global1\_part6.sh

```bash
chmod +x Scripts_bash/run_global1_part6.sh
./Scripts_bash/run_global1_part6.sh
```

### Step 7: run\_global1\_part456val.sh 

```bash
chmod +x Scripts_bash/run_global1_part456val.sh
./Scripts_bash/run_global1_part456val.sh
```

### Step 8: run\_global1\_part456test.sh 

```bash
chmod +x Scripts_bash/run_global1_part456test.sh
./Scripts_bash/run_global1_part456test.sh
```

---

## Data Availability

Please note that the dataset referenced here is **not located in the same directory** as mentioned above.
The dataset is available upon request.
Please contact:

```
emiliatugol@gmail.com
```
## PygoSemSim

Please refer to the following link for PygoSemSim:

https://github.com/mojaie/pygosemsim

## How to cite this work
Coming soon 
