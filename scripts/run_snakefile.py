# Time:        2023-08-02
# Author:      HaosenLuo
# FileName:    run_snakefile.py
# 制作样本信息的snakefile
from optparse import OptionParser
from pathlib import Path
import os


def run_snakefile(project, chip_id, analysis_id, sample_id, sample_name,
                  ref_genome_version, output_snakefile,
                  common_chrom_dup_threshold='2.8-3.2',
                  common_chrom_del_threshold='0.8-1.2',
                  male_y_chrom_dup_threshold='1.8-2.2',
                  male_y_chrom_del_threshold='0-0.4'):
    sample_snakefile = f"from pathlib import Path\n"\
                       f"from optparse import OptionParser\n"\
                       f"from datetime import datetime\n" \
                       f"import os\n"\
                       f"RESULT_PATH='/ifs1/result'\n"\
                       f"RESULT_PROJECT_PATH=Path(RESULT_PATH) / '{project}' / 'process'\n"\
                       f"if '{project}' == 'PGS_1M':\n"\
                       f"\tCHIP_ID_ANALYSIS_ID='_auto_'.join(['{chip_id}', '{analysis_id}'])\n"\
                       f"elif '{project}' == 'NIPT':\n"\
                       f"\tCHIP_ID_ANALYSIS_ID='_'.join(['{chip_id}', '{analysis_id}'])\n"\
                       f"elif '{project}' == 'DEMUX':\n"\
                       f"\tCHIP_ID_ANALYSIS_ID='_'.join(['ASA', '{chip_id}', '{analysis_id}'])\n"\
                       f"else:\n"\
                       f"\tCHIP_ID_ANALYSIS_ID='_'.join(['{chip_id}', '{analysis_id}'])\n"\
                       f"SAMPLE_ID_SAMPLE_NAME='_'.join(['{sample_id}', '{sample_name}'])\n"\
                       f"TARGET_FILE = '{sample_id}.bed'\n"\
                       f"SAMPLE_PATH = RESULT_PROJECT_PATH / CHIP_ID_ANALYSIS_ID / SAMPLE_ID_SAMPLE_NAME / TARGET_FILE\n"\
                       f"COMMON_CHROMOSOME_DUP_RANGE='{common_chrom_dup_threshold}'\n"\
                       f"COMMON_CHROMOSOME_DEL_RANGE='{common_chrom_del_threshold}'\n"\
                       f"MALE_Y_CHROMOSOME_DUP_RANGE='{male_y_chrom_dup_threshold}'\n"\
                       f"MALE_Y_CHROMOSOME_DEL_RANGE='{male_y_chrom_del_threshold}'\n"\
                       f"REFERENCE_GENOME_VERSION='{ref_genome_version}'\n" \
                       f"GENDER_FILE_PATH=RESULT_PROJECT_PATH / CHIP_ID_ANALYSIS_ID / SAMPLE_ID_SAMPLE_NAME / '*.gender'\n" \
                       # f"GENDER_GET=os.popen('head -1 {str(GENDER_FILE_PATH)}').read().strip()"\
                       # f"GENDER_OPTION=GENDER_GET\n"
    # Snakefile合并
    base_snakefile_path = '/ifs1/home/luohaosen/project/cnv_pathogenic_anno/Snakefile'
    Path(output_snakefile).parent.mkdir(exist_ok=True, parents=True)
    with open(output_snakefile, 'w') as out_writer:
        out_writer.write(sample_snakefile + '\n' + open(base_snakefile_path, 'r').read())
    output_cnv_anno_tsv = Path(output_snakefile).parent / f"{sample_id}.tsv"
    # 运行含样本信息的Snakefile
    os.system(f"snakemake -s {output_snakefile} --cores 1 {output_cnv_anno_tsv}")


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--project-opt', '--p', dest='project', help="The sample for project")
    parser.add_option('--chip_id', '--cid', dest='chip_id', help="The chip id of sample")
    parser.add_option('--analysis_id', '--aid', dest='analysis_id', help="The analysis id of sample")
    parser.add_option('--sample_id', '--sid', dest='sample_id', help="The sample id of sample")
    parser.add_option('--sample_name', '--nid', dest='sample_name', help="The sample name of sample")
    # parser.add_option('--gender', '--sex-gender', dest='gender', help="man or woman")
    parser.add_option('--ref_version', dest='ref_genome_version', help="The reference genome version like hg19 or hg38")
    parser.add_option('--output_snakefile', '--out', dest='output_snakefile', help="The output file of sample info")
    parser.add_option('--pc', '--dup-common-threshold', dest='common_chrom_dup_threshold',
                      help="Like dup threshold min-max", default='2.8-3.2')
    parser.add_option('--lc', '--del-common-threshold', dest='common_chrom_del_threshold',
                      help="Like del threshold min-max", default='0.8-1.2')
    parser.add_option('--pm', '--dup-man-threshold', dest='male_y_chrom_dup_threshold',
                      help="Like dup threshold min-max for man", default='1.8-2.2')
    parser.add_option('--lm', '--del-man-threshold', dest='male_y_chrom_del_threshold',
                      help="Like del threshold min-max for man", default='0.0-0.4')
    options, args = parser.parse_args()
    run_snakefile(project=options.project, chip_id=options.chip_id, analysis_id=options.analysis_id,
                  sample_id=options.sample_id, sample_name=options.sample_name,
                  # gender=options.gender,
                  ref_genome_version=options.ref_genome_version, output_snakefile=options.output_snakefile,
                  common_chrom_dup_threshold=options.common_chrom_dup_threshold,
                  common_chrom_del_threshold=options.common_chrom_del_threshold,
                  male_y_chrom_dup_threshold=options.male_y_chrom_dup_threshold,
                  male_y_chrom_del_threshold=options.male_y_chrom_del_threshold)
