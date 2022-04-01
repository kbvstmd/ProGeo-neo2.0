#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2021/10/9 18:48
# @Author  : Yu Zhang
import pandas as pd
import numpy as np
import re
import os


class NeoantigenUtils:

    def __init__(self, file_path1, file_path2=None):
        self.temp_file = file_path1
        self.out_file = file_path2

    @staticmethod
    def activate_env():
        os.system("source ~/.bashrc")

    @staticmethod
    def run_command(command, show=True):
        for c in command:
            if show:
                print(c)
            os.system(c)

    @staticmethod
    def create_dirs(dir_name):
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)

    @staticmethod
    def read_fasta():
        # 读取参考序列文件
        reuni_df0 = pd.read_csv("reference_file/reuniprot.fasta", sep='\t', header=None,
                                names=['gene', 'gene2', 'sequence'])
        reuni_df1 = reuni_df0[['gene', 'sequence']]
        reuni_df2 = reuni_df0[['gene2', 'sequence']]
        reuni_df2 = reuni_df2.rename(columns={'gene2': 'gene'})  # 修改列名
        reuni_df = pd.concat([reuni_df1, reuni_df2])  # 合并
        reuni_df.drop_duplicates(inplace=True)  # 去重
        reuni_df.reset_index(drop=True, inplace=True)  # 重置索引
        return reuni_df

    @staticmethod
    def sub_seq(seq_str, pos, pep_len):
        """

        :param seq_str: original sequence
        :param pos: mutation
        :param pep_len: peptide lengths
        :return: peptides list
        """
        pep_list = []
        for i in range(pep_len):
            pep_str = seq_str[pos - (pep_len - i):pos + i]
            if len(pep_str) == pep_len:
                mut_pos = pep_len - i
                pep_list.append([mut_pos, pep_str])
        return pep_list

    @staticmethod
    def find_gene_c_p(sample_df, pattern_c, pattern_p):
        result_list = []
        for sample in sample_df['AAChange.refGene']:
            str_list = sample.split(',')
            for s in str_list:
                gene = s.split(':')[0]
                c = re.findall(pattern_c, s)
                p = re.findall(pattern_p, s)
                if len(c) and len(p) and len(c) == len(p):
                    for i in range(len(c)):
                        result_list.append([gene, c[i], p[i]])
        result_df = pd.DataFrame(result_list, columns=['gene', 'cDNA', 'protein'])
        result_df.drop_duplicates(inplace=True)  # 去重
        result_df.reset_index(drop=True, inplace=True)  # 重置索引
        return result_df

    @staticmethod
    def alter_seq(seq, protein):  # 改变氨基酸字符
        alt_bef = protein[0]
        num = int(protein[1:-1])
        return '' if num > len(seq) or seq[num - 1] != alt_bef else seq[:num - 1] + protein[-1] + seq[num:]

    @staticmethod
    def alter_seq_del(seq, protein):
        if '_' in protein:
            str_list = protein.split('_')
            a1 = str_list[0][0]
            x1 = int(str_list[0][1:])
            x2 = int(str_list[1][1:])
            return '' if x2 > len(seq) or seq[x1 - 1] != a1 else seq[:x1 - 1] + seq[x2:]
        else:
            alt_bef = protein[0]
            num = int(protein[1:])
            return '' if num > len(seq) or seq[num - 1] != alt_bef else seq[:num - 1] + seq[num:]

    @staticmethod
    def alter_seq_ins(seq, protein):
        if 'delins' in protein:
            str_list = protein.split('delins')
            a1 = str_list[0][0]
            pos = int(str_list[0][1:])
            rep = str_list[1]
            return '' if pos > len(seq) or seq[pos - 1] != a1 else seq[:pos - 1] + rep + seq[pos:]
        else:
            str_list = protein.split('ins')
            ins = str_list[1]
            alt_bef = protein[0]
            pos = int(str_list[0].split('_')[0][1:])
            return '' if pos > len(seq) or seq[pos - 1] != alt_bef else seq[:pos] + ins + seq[pos:]

    @staticmethod
    def get_aac(protein, seq):
        pos = int(protein.split('_')[0][1:])
        if pos < len(seq):
            return protein[0] + seq[pos - 1]
        else:
            return protein[0]

    @staticmethod
    def get_aac_ins(protein, seq):
        if 'delins' in protein:
            str_list = protein.split('delins')
            return protein[0] + str_list[1][0]
        else:
            str_list = protein.split('_')[1].split('ins')
            return str_list[0][0] + str_list[1][0]

    def create_peps(self, alt_df, min_len, max_len):
        mhc_df = pd.DataFrame()
        for l in range(min_len, max_len + 1):
            mhc_list = []
            for i in range(alt_df.shape[0]):
                gene, c, p, mut, seq, alt_seq, pos, AAc = alt_df.iloc[i].values.tolist()
                pep_list = self.sub_seq(alt_seq, pos, l)
                for mut_pos, pep in pep_list:
                    mhc_list.append([gene, c, p, AAc, mut_pos, pep])
            temp = pd.DataFrame(mhc_list, columns=['gene', 'cDNA', 'protein', 'AAc', 'mut_pos', 'peptide'])
            mhc_df = pd.concat([mhc_df, temp])
        return mhc_df

    @staticmethod
    def output_peps(STR_MHC_df: pd.DataFrame):
        """
        :param STR_MHC_df: NM & peptide pairs
        :return: only one column with NM or peptide
        """
        STR_MHC_df.drop_duplicates(inplace=True)  # 去重
        # STR_MHC_df.reset_index(drop=True, inplace=True)  # 重置索引
        STR_MHC_list = STR_MHC_df.values.tolist()
        STR_MHC_list = np.concatenate(STR_MHC_list)
        final_MHC_df = pd.DataFrame(STR_MHC_list, columns=['peptides'])
        return final_MHC_df

    def save_fasta(self, mhc_df, tag, save=1):
        mhc_df['title'] = mhc_df.apply(lambda row: '>' + row['gene'] + ':' + row['AAc'] + ':' + str(row['mut_pos']),
                                       axis=1)
        if save == 1:  # 保存Ⅰ类MHC
            STR_MHC_df = self.output_peps(mhc_df[['title', 'peptide']].copy())
            STR_MHC_df.to_csv('outfile-wes/fasta_files/' + tag + '_MHC_1.fasta', index=False,
                              header=None, sep='\t')
        else:  # 保存Ⅱ类MHC
            for length in range(15, 31):
                STR_MHC_df = mhc_df[mhc_df['peptide'].str.len() == length]
                STR_MHC_df = self.output_peps(STR_MHC_df[['title', 'peptide']].copy())
                # 存入fasta文件
                STR_MHC_df.to_csv('outfile-wes/fasta_files/' + tag + '_MHC_2_' + str(length) + '.fasta',
                                  index=False, header=None, sep='\t')

    @staticmethod
    def save_csv(mhc_df, tag, save=1):
        if save == 1:  # 保存Ⅰ类MHC
            mhc_df[['gene', 'AAc', 'mut_pos', 'peptide', 'cDNA', 'protein']].to_csv(
                'outfile-wes/csv_files/' + tag + '_MHC_1.csv', index=False)
        else:  # 保存Ⅱ类MHC
            for length in range(15, 31):
                STR_MHC_df = mhc_df[mhc_df['peptide'].str.len() == length]
                # 存入fasta文件
                STR_MHC_df[['gene', 'AAc', 'mut_pos', 'peptide', 'cDNA', 'protein']].to_csv(
                    'outfile-wes/csv_files/' + tag + '_MHC_2_' + str(length) + '.csv', index=False)

    def handle1(self, sample, reuniprot, sub_pep='y'):  # SNV
        tag = 'snv'
        sample_snv = sample[sample['ExonicFunc.refGene'].isin(['nonsynonymous SNV'])]
        sample_snv.reset_index(drop=True, inplace=True)
        pattern_c = re.compile(r'c.([A-Z]\d{1,8}[A-Z]):')
        pattern_p = re.compile(r'p.([A-Z]\d{1,8}[A-Z])')
        result_df = self.find_gene_c_p(sample_snv, pattern_c, pattern_p)
        result_df['mut'] = 'nonsynonymous SNV'
        # 保存为txt文件
        result_df.to_csv(self.out_file + '/txt_files/nonsynonymous_SNV.txt', sep='\t', index=False)
        # 根据gene名查找对应序列
        merge_snv_df = pd.merge(result_df, reuniprot, how='inner', left_on='gene', right_on='gene')
        # 在突变位点替换氨基酸
        alt_snv_df = merge_snv_df.copy()
        alt_snv_df['alt_seq'] = alt_snv_df[['sequence', 'protein']].apply(
            lambda row: self.alter_seq(row['sequence'], row['protein']), axis=1)
        alt_snv_df = alt_snv_df[~alt_snv_df['alt_seq'].isin([''])]  # 删除空al_seq的行
        alt_snv_df['pos'] = alt_snv_df['protein'].apply(lambda p: int(p[1:-1]))
        alt_snv_df['AAc'] = alt_snv_df['protein'].apply(lambda p: p[0] + p[-1])
        alt_snv_df.to_csv(self.out_file + '/csv_files/alt_' + tag + '.csv', index=False)
        if sub_pep == 'y':
            self.peps_save(alt_snv_df, 'snv')
            print("snv MHC files have been saved.")

    def handle2(self, sample, reuniprot, sub_pep='y'):  # DEL,缺失
        tag = 'del'
        sample_del = sample[sample['ExonicFunc.refGene'].isin(['nonframeshift deletion'])]
        sample_del.reset_index(drop=True, inplace=True)
        pattern_c = re.compile(r'c.(\d{1,4}_\d{1,4})del:')
        pattern_p = re.compile(r'p.([A-Z].*?)del')
        result_df = self.find_gene_c_p(sample_del, pattern_c, pattern_p)
        result_df['mut'] = 'nonframeshift deletion'
        # 保存为txt文件
        result_df.to_csv(self.out_file + '/txt_files/nonframeshift_deletion.txt', sep='\t', index=False)
        # 根据gene名查找对应序列
        merge_del_df = pd.merge(result_df, reuniprot, how='inner', left_on='gene', right_on='gene')
        # 在突变位点替换氨基酸
        alt_del_df = merge_del_df.copy()
        alt_del_df['alt_seq'] = alt_del_df[['sequence', 'protein']].apply(
            lambda row: self.alter_seq_del(row['sequence'], row['protein']), axis=1)
        alt_del_df = alt_del_df[~alt_del_df['alt_seq'].isin([''])]
        alt_del_df['pos'] = alt_del_df['protein'].apply(lambda p: int(p.split('_')[0][1:]))
        alt_del_df['AAc'] = alt_del_df.apply(lambda row: self.get_aac(row['protein'], row['alt_seq']), axis=1)
        alt_del_df.to_csv(self.out_file + '/csv_files/alt_' + tag + '.csv', index=False)
        if sub_pep == 'y':
            self.peps_save(alt_del_df, 'del')
            print("deletion MHC files have been saved.")

    def handle3(self, sample, reuniprot, sub_pep='y'):  # 缺失插入, delins
        tag = 'ins'
        sample_delins = sample[sample['ExonicFunc.refGene'].isin(['nonframeshift insertion'])]  # InDel=DEL&INS
        sample_delins.reset_index(drop=True, inplace=True)
        pattern_c = re.compile(r'c.(\d{1,4}_\d{1,4}ins[A-Z]{1,}):')
        pattern_p = re.compile(r'p.(.*[A-Z]{1,}$)')
        result_df = self.find_gene_c_p(sample_delins, pattern_c, pattern_p)
        result_df['mut'] = 'nonframeshift insertion'
        # 保存为txt文件
        result_df.to_csv(self.out_file + '/txt_files/nonframeshift_insertion.txt', sep='\t', index=False)
        # 根据gene名查找对应序列
        merge_delins_df = pd.merge(result_df, reuniprot, how='inner', left_on='gene', right_on='gene')
        # 在突变位点 插入/替换 氨基酸
        alt_ins_df = merge_delins_df.copy()
        alt_ins_df['alt_seq'] = alt_ins_df[['sequence', 'protein']].apply(
            lambda row: self.alter_seq_ins(row['sequence'], row['protein']), axis=1)
        alt_ins_df = alt_ins_df[~alt_ins_df['alt_seq'].isin([''])]
        alt_ins_df['pos'] = alt_ins_df['protein'].apply(
            lambda p: int(p.split('delins')[0][1:]) if 'delins' in p else int(p.split('_')[1].split('ins')[0][1:]))
        alt_ins_df['AAc'] = alt_ins_df.apply(lambda row: self.get_aac_ins(row['protein'], row['alt_seq']), axis=1)
        alt_ins_df.to_csv(self.out_file + '/csv_files/alt_' + tag + '.csv', index=False)
        if sub_pep == 'y':
            self.peps_save(alt_ins_df, 'ins')
            print("insert MHC files have been saved.")

    def handle4(self, sample, reuniprot, sub_pep='y'):  # 移码突变, fs
        tag = 'fs'
        sample_fs_del = sample[sample['ExonicFunc.refGene'].isin(['frameshift deletion'])]
        sample_fs_ins = sample[sample['ExonicFunc.refGene'].isin(['frameshift insertion'])]
        pattern_c = re.compile(r'c.(\d{1,4}_\d{1,4})del:')
        pattern_c2 = re.compile(r'c.(\d{1,4}_\d{1,4}ins[A-Z]{1,}):')
        pattern_p = re.compile(r'p.([A-Z]\d{1,4}[A-Z])fs*')
        result_df1 = self.find_gene_c_p(sample_fs_del, pattern_c, pattern_p)
        result_df1['mut'] = 'frameshift deletion'
        result_df2 = self.find_gene_c_p(sample_fs_ins, pattern_c2, pattern_p)
        result_df2['mut'] = 'frameshift insertion'
        result_df = pd.concat([result_df1, result_df2])
        result_df.reset_index(drop=True, inplace=True)
        # 保存为txt文件
        result_df.to_csv(self.out_file + '/txt_files/frameshift_fusion.txt', sep='\t', index=False)
        # 根据gene名查找对应序列
        merge_fs_df = pd.merge(result_df1, reuniprot, how='inner', left_on='gene', right_on='gene')
        # 在突变位点替换氨基酸信息
        alt_fs_df = merge_fs_df.copy()
        alt_fs_df['alt_seq'] = alt_fs_df[['sequence', 'protein']].apply(
            lambda row: self.alter_seq(row['sequence'], row['protein']), axis=1)
        alt_fs_df = alt_fs_df[~alt_fs_df['alt_seq'].isin([''])]  # 删除空al_seq的行
        alt_fs_df['pos'] = alt_fs_df['protein'].apply(lambda p: int(p[1:-1]))
        alt_fs_df['AAc'] = alt_fs_df['protein'].apply(lambda p: p[0] + p[-1])
        alt_fs_df.to_csv(self.out_file + '/csv_files/alt_' + tag + '.csv', index=False)
        if sub_pep == 'y':
            self.peps_save(alt_fs_df, 'fs')
            print("fs MHC files have been saved.")

    def save_var_pro(self):  # 保存突变的长肽序列
        alt_df = pd.DataFrame()
        mut_types = ['snv', 'del', 'ins', 'fs']
        for mut_type in mut_types:
            mut_df = pd.read_csv(self.out_file + '/csv_files/alt_' + mut_type + '.csv')
            alt_df = pd.concat([alt_df, mut_df])
        alt_df['title'] = alt_df.apply(lambda row: '>' + row['gene'] + '|p.' + row['protein'], axis=1)
        var_pro_df = self.output_peps(alt_df[['title', 'alt_seq']].copy())
        var_pro_df.to_csv(self.out_file + '/Varsequence.fasta', index=False, header=None, sep='\t')
        print("Varsequence.fasta files have been saved.")

    # def sub_21_seq(self, alt_seq, pos):
    #     pos = pos - 1
    #     if pos - 10 < 0:
    #         return alt_seq[:pos + 11 + (10 - pos)]
    #     elif pos + 11 > len(alt_seq):
    #         return alt_seq[-21:]
    #     return alt_seq[pos - 10:pos + 11]
    #
    # def save_var_pro_short(self, index):  # 保存21-mer短肽序列
    #     alt_df = pd.DataFrame()
    #     mut_types = ['snv', 'del', 'ins', 'fs']
    #     for mut_type in mut_types:
    #         mut_df = pd.read_csv('outputs/sample' + str(index) + '/csv_files/alt_' + mut_type + '.csv')
    #         alt_df = pd.concat([alt_df, mut_df])
    #     alt_df['title'] = alt_df.apply(lambda row: '>' + row['gene'] + '|p.' + row['protein'], axis=1)
    #     alt_df['short_seq'] = alt_df.apply(lambda row: self.sub_21_seq(row['alt_seq'], row['pos']), axis=1)
    #     var_pro_df = self.output_peps(alt_df[['title', 'short_seq']])
    #     var_pro_df.to_csv('outputs/sample' + str(index) + '/fasta_files/VarproSeq_s.fasta',
    #                       index=False, header=None, sep='\t')

    def peps_save(self, alt_df, tag):
        # 保存短肽为fasta文件和csv文件
        mhc_df1 = self.create_peps(alt_df, 8, 11)
        self.save_fasta(mhc_df1, tag)
        self.save_csv(mhc_df1, tag)
        mhc_df2 = self.create_peps(alt_df, 15, 30)
        self.save_fasta(mhc_df2, tag, 2)
        self.save_csv(mhc_df2, tag, 2)
