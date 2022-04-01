import os
import pandas as pd
import re
from neoantigen_utils import NeoantigenUtils


def neo_pred_I():
    # input the hla allele:  HLA-A03:01,HLA-B07:02,HLA-B35:03,HLA-C07:02,HLA-C04:01
    hla_allele1 = input(
        "please input an HLA class I allele like 'HLA-A03:01' or multiple alleles like 'HLA-A03:01,HLA-B07:02,HLA-B35:03':")
    hla_allele1 = hla_allele1.replace(' ', '')
    tags = mut_types
    ### MHC I:
    command = []
    for tag in tags[:-1]:
        command.append(
            "netMHCpan -f outfile-wes/fasta_files/" + tag + "_MHC_1.fasta -a " + hla_allele1 + " -s -BA > " + utils.temp_file + "/MHC-I/" + tag + "/" + tag + "_1.out &")
    # fusion:
    command.append(
        "netMHCpan -f outfile-rna/fusion/fusion-pep.fasta -a " + hla_allele1 + " -s -BA > " + utils.temp_file + "/MHC-I/fusion/fusion_1.out &")
    utils.run_command(command, False)


def neo_pred_II():
    # input the hla allele:
    hla_allele2 = input(
        "please input an HLA class II allele like 'DRB1_0101' or multiple alleles like 'DRB1_0102,DRB1_0301,DQB1_0101':")
    hla_allele2 = hla_allele2.replace(' ', '')
    tags = mut_types[:-1]
    ### MHC II:
    command = []
    for tag in tags:
        for i in range(15, 31):
            command.append(
                "netMHCIIpan -f outfile-wes/fasta_files/" + tag + "_MHC_2_" + str(i) + ".fasta -a " + hla_allele2 + " -BA -s -length " + str(i) + " > " + utils.temp_file + "/MHC-II/" + tag + "/" + tag + "_2_" + str(i) + "_fasta.out &")
    utils.run_command(command, False)


def mutant_judge(peptide, str0, split_ch):
    str_list = str0.split(split_ch)
    if len(str_list) == 3 and len(str_list[1]) == 2:
        aachange = str_list[1]
        pos = int(str_list[2])
        if pos <= len(peptide):
            return 1 if peptide[pos - 1] == aachange[1] else 0
        return 0
    return 0


def netmhc_filter(input1, hla_class=1):
    """
    过滤netMHCpan的结果
    :param input1: netMHCpan的输出文件
    :param hla_class: HLA类型
    :return:
    """
    col_names = ''
    split_ch = '_' if hla_class == 1 else ':'
    for line in input1:  # 找出第一个带有行名的位置索引,作为列名
        if 'Pos' in line:
            col_names = re.split(r"[ ]+", line)
            break
    out_data = []
    for line in input1:  # 将所有信息放进list以便转换df
        s_list = re.split(r"[ ]+", line)
        if 'WB\n' in s_list[-1] or 'SB\n' in s_list[-1]:
            out_data.append(s_list)
    input1.close()  # 关闭文件
    if hla_class == 1:
        col_names.insert(-1, 'Bind')
    # 生成dataframe
    out_df = pd.DataFrame(out_data, columns=col_names)
    out_df = out_df.dropna(how='any')  # 删除有空值的行
    out_df['BindLevel'] = out_df['BindLevel\n'].apply(lambda bind: 'WB' if 'WB' in bind else 'SB')  # 提取BindLevel
    if hla_class == 2:
        out_df.rename(columns={'Affinity(nM)': 'Aff(nM)'}, inplace=True)  # 修改列名
    identity_list = out_df['Identity'].tolist()
    if identity_list == []:
        return None
    elif '--' in identity_list[0]:  # 融合基因
        output1 = out_df.copy()
        output1 = output1.rename(columns={
            'MHC': 'HLA',
            'Identity': 'gene'
        })
        select_col = ['HLA', 'gene', 'Peptide', 'Aff(nM)', '%Rank_BA', 'BindLevel']
    else:  # 判断突变信息
        out_df['judge'] = out_df.apply(lambda row: mutant_judge(row['Peptide'], row['Identity'], split_ch), axis=1)
        output1 = out_df[out_df['judge'] == 1].copy()
        output1['gene'] = output1['Identity'].apply(lambda x: x.split(split_ch)[0])
        output1['AAchange'] = output1['Identity'].apply(lambda x: x.split(split_ch)[1])
        output1['mut_pos'] = output1['Identity'].apply(lambda x: int(x.split(split_ch)[2]))
        output1['judge_pos'] = output1[['Pos', 'mut_pos']].apply(
            lambda row: 0 if int(row['Pos']) > row['mut_pos'] else 1, axis=1)
        output1 = output1[output1['judge_pos'] == 1]
        output1.rename(columns={'MHC': 'HLA'}, inplace=True)
        select_col = ['HLA', 'gene', 'Peptide', 'AAchange', 'mut_pos', 'Aff(nM)', '%Rank_BA', 'BindLevel']
    out = output1[select_col].copy()
    out.drop_duplicates(inplace=True)
    return out


def neo_candid(hla_class=1):
    """
    处理netMHCpan输出的文件，得到候选新抗原
    I类：直接处理后保存
    II类：15~30-mer肽段分开处理，合并为一个txt文件后保存
    :param:
    :return:
    """
    mut_type = mut_types if hla_class == 1 else mut_types[:-1]
    if hla_class == 1:
        for tag in mut_type:
            input1 = open(utils.temp_file + "/MHC-I/" + tag + "/" + tag + "_1.out", "r")  # 读取文件
            out = netmhc_filter(input1)
            if out is not None:
                out.to_csv(utils.out_file + "/MHC-I/" + tag + "/I-" + tag + "-candidate-neoantigens.txt", index=False, sep='\t')
    else:
        for tag in mut_type:
            final_out = pd.DataFrame()
            for i in range(15, 31):
                input2 = open(utils.temp_file + "/MHC-II/" + tag + "/" + tag + "_2_" + str(i) + "_fasta.out", "r")  # 读取II类文件
                out = netmhc_filter(input2, 2)
                if out is not None:
                    final_out = pd.concat([final_out, out])
            final_out.to_csv(utils.out_file + "/MHC-II/" + tag + "/II-" + tag + "-candidate-neoantigens.txt", index=False, sep='\t')


if __name__ == '__main__':
    utils = NeoantigenUtils("outfile3", "outfile-candidate-neoantigens")
    directs = ['MHC-I', 'MHC-II']
    mut_types = ['snv', 'del', 'ins', 'fs', 'fusion']
    utils.create_dirs(utils.temp_file)
    utils.create_dirs(utils.out_file)
    for d in directs:
        utils.create_dirs(utils.temp_file + '/' + d)
        directs2 = mut_types if d == 'MHC-I' else mut_types[:-1]
        for d2 in directs2:
            utils.create_dirs(utils.temp_file + '/' + d + '/' + d2)
    for d in directs:
        utils.create_dirs(utils.out_file + '/' + d)
        directs2 = mut_types if d == 'MHC-I' else mut_types[:-1]
        for d2 in directs2:
            utils.create_dirs(utils.out_file + '/' + d + '/' + d2)
    neo_pred_I()
    neo_pred_II()
    neo_candid()
    neo_candid(2)
