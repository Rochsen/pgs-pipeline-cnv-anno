# Time:        2023-08-03
# Author:      HaosenLuo
# FileName:    pgs_cnv_anno.py
import os
import re
from pathlib import Path
from optparse import OptionParser
import time
import pandas as pd

CNV_TO_BED_PY = '/ifs1/home/luohaosen/project/cnv_pathogenic_anno/scripts/cnv_to_bed.py'
ANNO_CLASSIFY_CNV_PY = '/ifs1/home/luohaosen/project/cnv_pathogenic_anno/scripts/anno_cnv_classify.py'


# 获取某个批次的汇总文件
def __get_run_sample_file(chip_id, analysis_id, project):
    if project == 'PGS_1M':
        result_merge_file = '_'.join([chip_id, 'result_merge_new.tsv'])
        return Path('/ifs1/result/') / project / 'process' / '_auto_'.join([chip_id, analysis_id]) / result_merge_file


# 主要流程（注意判断是否存在bed文件，不存在就跳过）
def run_cnv_anno(chip_id, analysis_id,
                 project='PGS_1M', ref_genome_version='hg19',
                 common_chrom_dup_threshold='2.8', common_chrom_del_threshold='1.2',
                 male_y_chrom_dup_threshold='1.8', male_y_chrom_del_threshold='0.4',
                 ):
    input_chip_id = chip_id
    input_analysis_id = analysis_id
    input_project = project
    ref_version = ref_genome_version
    run_file = __get_run_sample_file(input_chip_id, input_analysis_id, input_project)
    data = pd.read_table(run_file)
    for _, row in data.iterrows():
        sample_id = row['Sample_id']
        sample_name = row['Sample_name']
        primer_label = row['引物标签']
        sample_parent_dir = Path('/ifs1/result/') / project / 'process' / '_auto_'.join([chip_id, analysis_id])
        sample_directory = '_'.join([sample_id, sample_name])
        sample_input_bed = sample_parent_dir / sample_directory / f"{sample_id}.bed"
        sample_output_bed = sample_parent_dir / sample_directory / f"{sample_id}_run_cnv/{sample_id}.bed"
        sample_output_tsv = sample_parent_dir / sample_directory / f"{sample_id}_run_cnv/{sample_id}.tsv"
        # sample_classify_cnv_result_path = sample_parent_dir / sample_directory / "run_cnv/classifycnv_anno/"
        # 检查bed文件是否存在
        if sample_input_bed.exists():
            # 查找性别
            gender_file = sample_parent_dir / sample_directory / f"{sample_id}_{input_chip_id}_L01_{primer_label}.gender"
            if gender_file.exists():
                gender_str = os.popen(f"head -1 {gender_file}").read().strip()
                # 构建适合于ClassifyCnv.py的bed文件
                os.system(f"python {CNV_TO_BED_PY} --i {sample_input_bed}\
                                                   --o {sample_output_bed}\
                                                   --gender {gender_str}\
                                                   --dup-common-threshold {common_chrom_dup_threshold}\
                                                   --del-common-threshold {common_chrom_del_threshold}\
                                                   --dup-man-threshold {male_y_chrom_dup_threshold}\
                                                   --del-man-threshold {male_y_chrom_del_threshold}")
                # 跑ClassifyCNV.py，获取注释信息
                if os.popen(f"less {sample_output_bed}").read().strip() != '':
                    os.system(f"python {ANNO_CLASSIFY_CNV_PY} --i {sample_output_bed}\
                                                              --o {sample_output_tsv}\
                                                              --g {ref_version}")
                else:
                    print(f"{sample_output_tsv} is not Exists")
            else:
                print(f"{gender_file} not exists")
        else:
            print(f"{sample_input_bed} not exists")


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--chip_id', '--cid', dest='chip_id', help="The chip id of sample")
    parser.add_option('--analysis_id', '--aid', dest='analysis_id', help="The analysis id of sample")
    parser.add_option('--project-opt', dest='project', default='PGS_1M', help="The sample for project")
    parser.add_option('--ref_version', dest='ref_genome_version', default='hg19',
                      help="The reference genome version like hg19 or hg38")
    parser.add_option('--pc', '--dup-common-threshold', dest='common_chrom_dup_threshold',
                      help="Like dup threshold min-max", default='2.8')
    parser.add_option('--lc', '--del-common-threshold', dest='common_chrom_del_threshold',
                      help="Like del threshold min-max", default='1.2')
    parser.add_option('--pm', '--dup-man-threshold', dest='male_y_chrom_dup_threshold',
                      help="Like dup threshold min-max for man", default='1.8')
    parser.add_option('--lm', '--del-man-threshold', dest='male_y_chrom_del_threshold',
                      help="Like del threshold min-max for man", default='0.4')
    options, args = parser.parse_args()
    run_cnv_anno(chip_id=options.chip_id,
                 analysis_id=options.analysis_id,
                 project=options.project,
                 ref_genome_version=options.ref_genome_version,
                 common_chrom_dup_threshold=options.common_chrom_dup_threshold,
                 common_chrom_del_threshold=options.common_chrom_del_threshold,
                 male_y_chrom_dup_threshold=options.male_y_chrom_dup_threshold,
                 male_y_chrom_del_threshold=options.male_y_chrom_del_threshold,
                 )
