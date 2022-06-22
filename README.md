# ProGeo-neo v2.0

A One-stop software for Neoantigen Prediction and Filtering based on the Proteogenomics strategy.
## Neoantigen Prediction and Filtering Pipeline
![](pipeline.jpg)
## 1. Introduction
ProGeo-neo v2.0 is a one-stop software solution for predicting neoantigens from the paired tumour-normal WGS/WES data in FASTQ format. ProGeo-neo v2.0 provides new features such as in-frame indels, frameshift mutations, and gene fusion analysis. In addition, the new version supports the prediction of MHC class II-restricted neoantigens, i.e. peptides up to 30mer in length. The source of neoantigens has been expanded, allowing more candidate neoantigens to be identified for follow-up studies. In addition, we propose two more efficient screening approaches, mainly using an in-group authentic neoantigen database for neoantigen screening. The range of candidate peptides was effectively narrowed down to those that are more likely to elicit an immune response, providing a more meaningful reference for subsequent experimental validation. Compared to ProGeo-neo, the ProGeo-neo v2.0 performs well on the same data set, and this means that our upgrades are necessary, both in updated functionality and improved accuracy. 
## 2. Running environment
ProGeo-neo2.0 requires a Linux operation system(centos7) with Python(V3.7), Perl(V5.26) and Java(V1.7) installed.
## 3. Usage
Run the following codes before getting started.
```
cd ProGeo-neo2.0
bash start.sh 
# Users with root privileges can ignore the following:
chmod 755 software/bwa/bwa 
chmod 755 software/samtools/samtools 
chmod 755 software/bcftools/bcftools
```
### 3.1. WGS/WES processing
```
python wes_mutation_peptides.py /path/to/wes-tumor /path/to/wes-normal
# eg:
python wes_mutation_peptides.py test/wes/tumor test/wes/normal
```
Temporary files and final result files generated during data processing will be placed under outfile1 and outfile-wes directories respectively.
### 3.2. RNA processing
```
python RNA_seq_mutation.py /path/to/rna
# eg:
python RNA_seq_mutation.py test/rna
```
**notes:**  
HLA types are required for predicting neoantigens. You can either enter your own types or infer HLA types by using the following commands.
Input "y" if users need HLA types from RNA sequences when the system prompts: "Predicting HLA class I types from next-generation sequencing data: (y/n)?" or "Predicting HLA class II types from next-generation sequencing data: (y/n)?", otherwise input "n".

Temporary files and final result files generated during data processing will be placed under the outfile2 and the outfile-rna respectively. It should be noted that the predicted HLA types will be stored in .tsv files under the outfile-rna/hla directory, contributing to subsequent neoantigen screening strategies. In addition, the synthetic mutant long peptides which will be saved in the Varsequence.fasta will be used for mass spectrometry library construction later.
### 3.3. Neoantigen prediction
```
python neoantigen_prediction.py
```
**notes:**  
Input the HLA types predicted in 3.2 or other types that the user interested in when the system prompts: 
"please input an HLA class I allele like 'HLA-A03:01' or multiple alleles like 'HLA-A03:01,HLA-B07:02,HLA-B35:03':" and "please input an HLA class II allele like 'DRB1_0101' or multiple alleles like 'DRB1_0102,DRB1_0301,DQB1_0101':".

The predicted HLA class I and HLA class II candidate neoantigens will be stored in the MHC-I and MHC-II catalogs under the outfile-candidate-neoantigens directory, with the temporary files under the outfile3 directory.
### 3.4. Mass spectrometry filtration
```
python MS_filtration.py /path/to/mass
# eg:
python MS_filtration.py test/mass
```
The result peptides filtered by this step will be stored in the filter-mass/Maxpep.txt file and be used for further neoantigen filtration.
### 3.5. Neoantigen filtration
```
python neoantigen_filtration.py
```
After different levels of strict threshold filtration, the final candidate neoantigens will be saved in the MHC-I and MHC-II directories under the outfile-filter-neoantigens directory.

For detailed software installation and usage of ProGeo-neo2.0, please read the User's Manual.

### 3.6. Cite
Liu, C.; Zhang, Y.; Jian, X.; Tan, X.; Lu, M.; Ouyang, J.; Liu, Z.; Li, Y.; Xu, L.; Chen, L.; Lin, Y.; Xie, L. ProGeo-Neo v2.0: A One-Stop Software for Neoantigen Prediction and Filtering Based on the Proteogenomics Strategy. Genes 2022, 13, 783.
