#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2021/12/24 16:01
# @Author  : Yu Zhang
"""
wes data processing(tumor/normal):
"""
import os
import sys
import pandas as pd
from neoantigen_utils import NeoantigenUtils


def handle01_bwa():
    command = ['bwa index reference_file/hg38/hg38.fa']
    utils.run_command(command)


def handle02_sam(dir_tumor, dir_normal):
    """
    Map reads to reference genome  BWA（sam）
    """
    tumor_files = os.listdir(dir_tumor)
    normal_files = os.listdir(dir_normal)
    if len(tumor_files) == 2 and len(normal_files) == 2:
        command = [
            r"bwa mem -t 4 -M -q -5 -R '@RG\tID:foo\tSM:bar\tLB:library1' reference_file/hg38/hg38.fa " + dir_tumor + "/" +
            tumor_files[0] + " " + dir_tumor + "/" + tumor_files[1] + " > " + utils.temp_file + "/Tumor-sample.sam",
            r"bwa mem -t 24 -M -q -5 -R '@RG\tID:foo\tSM:bar\tLB:library1' reference_file/hg38/hg38.fa " + dir_normal + "/" +
            normal_files[0] + " " + dir_normal + "/" + normal_files[1] + " > " + utils.temp_file + "/Adjacengt-sample.sam"]
        utils.run_command(command)
    else:
        print("please input files into the right path!")


def handle03_bam_sort():
    utils.activate_env()
    command = ["samtools fixmate -O bam " + utils.temp_file + "/Tumor-sample.sam " + utils.temp_file + "/Tumor-sample.bam",
               "samtools fixmate -O bam " + utils.temp_file + "/Adjacengt-sample.sam " + utils.temp_file + "/Adjacengt-sample.bam",
               "samtools sort -O bam -o " + utils.temp_file + "/Tumor_sorted_sample.bam -T temp " + utils.temp_file + "/Tumor-sample.bam",
               "samtools sort -O bam -o " + utils.temp_file + "/Adjacengt_sorted_sample.bam -T temp " + utils.temp_file + "/Adjacengt-sample.bam"]
    utils.run_command(command)


def handle04_mark_duplicates():
    command = [
        "gatk MarkDuplicates -I " + utils.temp_file + "/Tumor_sorted_sample.bam -O " + utils.temp_file + "/Tumor_sorted_markedup.bam -M " + utils.temp_file + "/Tumor_sorted_markedup_metric.bam",
        "gatk MarkDuplicates -I " + utils.temp_file + "/Adjacengt_sorted_sample.bam -O " + utils.temp_file + "/Adjacengt_sorted_markedup.bam -M " + utils.temp_file + "/Adjacengt_sorted_markedup_metric.bam"]
    utils.run_command(command)


def handle05_add_head():
    command = [
        "gatk AddOrReplaceReadGroups -I  " + utils.temp_file + "/Tumor_sorted_markedup.bam -O " + utils.temp_file + "/Tumor_sorted_markedup_Add.bam -ID 4 -LB lib1 -PL illumina -PU unit1 -SM cancer",
        "gatk AddOrReplaceReadGroups -I " + utils.temp_file + "/Adjacengt_sorted_markedup.bam -O " + utils.temp_file + "/Adjacengt_sorted_markedup_Add.bam -ID 4 -LB lib1 -PL illumina -PU unit1 -SM normal"]
    utils.run_command(command)


def handle06_recalibrate_bases():
    command = [
        "gatk BaseRecalibrator -I " + utils.temp_file + "/Tumor_sorted_markedup_Add.bam  -R  " + current_path + "/reference_file/hg38/hg38.fa  --known-sites reference_file/dbsnp_146.hg38.vcf -O " + utils.temp_file + "/Tumor_recal_data-sample.table",
        "gatk BaseRecalibrator -I " + utils.temp_file + "/Adjacengt_sorted_markedup_Add.bam  -R  " + current_path + "/reference_file/hg38/hg38.fa  --known-sites reference_file/dbsnp_146.hg38.vcf -O " + utils.temp_file + "/Adjacengt_recal_data-sample.table"]
    utils.run_command(command)


def handle07_ApplyBQSR():
    command = [
        "gatk ApplyBQSR -R reference_file/hg38/hg38.fa -I " + utils.temp_file + "/Tumor_sorted_markedup_Add.bam  --bqsr-recal-file  " + utils.temp_file + "/Tumor_recal_data-sample.table -O  " + utils.temp_file + "/Tumor_recal-sample.bam",
        "gatk ApplyBQSR -R reference_file/hg38/hg38.fa -I " + utils.temp_file + "/Adjacengt_sorted_markedup_Add.bam  --bqsr-recal-file  " + utils.temp_file + "/Adjacengt_recal_data-sample.table -O  " + utils.temp_file + "/Adjacengt_recal-sample.bam"]
    utils.run_command(command)


def handle08_somatic():
    command = [
        "nohup gatk Mutect2 -R reference_file/hg38/hg38.fa -I " + utils.temp_file + "/Adjacengt_recal-sample.bam -O " + utils.temp_file + "/Adjacengt_normal.vcf &",
        "nohup gatk Mutect2 -R reference_file/hg38/hg38.fa -I " + utils.temp_file + "/Tumor_recal-sample.bam -I " + utils.temp_file + "/Adjacengt_recal-sample.bam -tumor cancer -normal normal -pon " + utils.temp_file + "/Adjacengt_normal.vcf -O " + utils.temp_file + "/1_somatic_m2.vcf.gz -bamout " + utils.temp_file + "/2_tumor_normal_m2.bam &"]
    utils.run_command(command)


def handle09_gatk():
    """
    gatk FilterMutectCalls+vcftools（vcf）
    """
    command = [
        "gatk FilterMutectCalls -V " + utils.temp_file + "/1_somatic_m2.vcf.gz -R reference_file/hg38/hg38.fa -O " + utils.temp_file + "/somatic_m2.Filtered.vcf"]
    utils.run_command(command)


def handle10_annovar():
    command = [
        "perl software/annovar/table_annovar.pl " + utils.temp_file + "/somatic_m2.Filtered.vcf software/annovar/humandb1/ -buildver hg38 -out " + utils.out_file + "/wes_annovar_out -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput"]
    utils.run_command(command)


def handle11(sample, reuniprot):
    """
    create mutant peptides, save to csv./txt./fasta files
    """
    utils.handle1(sample, reuniprot)
    utils.handle2(sample, reuniprot)
    utils.handle3(sample, reuniprot)
    utils.handle4(sample, reuniprot)


if __name__ == '__main__':
    sys_args = sys.argv
    utils = NeoantigenUtils("outfile1", "outfile-wes")
    # current_path = os.path.abspath(__file__)
    current_path = os.path.abspath(os.path.dirname(__file__))
    # print(current_path)
    utils.create_dirs(utils.temp_file)
    utils.create_dirs(utils.out_file)
    dir_list = ['/txt_files', '/fasta_files', '/csv_files']
    for d in dir_list:
        utils.create_dirs(utils.out_file + d)
    utils.activate_env()
    handle01_bwa()
    handle02_sam(sys_args[1], sys_args[2])
    handle03_bam_sort()
    handle04_mark_duplicates()
    handle05_add_head()
    handle06_recalibrate_bases()
    handle07_ApplyBQSR()
    handle08_somatic()
    handle09_gatk()
    handle10_annovar()
    reuniprot = utils.read_fasta()
    # # 读取样本文件
    sample = pd.read_table(utils.out_file + '/wes_annovar_out.hg38_multianno.txt', sep="\t")
    sample = sample[sample["Otherinfo10"].str.contains("PASS")]
    handle11(sample, reuniprot)
