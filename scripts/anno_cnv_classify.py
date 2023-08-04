# Time:        2023-07-31
# Author:      HaosenLuo
# FileName:    anno_cnv_classify.py
import pandas as pd
from pathlib import Path
import os
from optparse import OptionParser
import sys
classify_cnv_py = '/ifs1/home/luohaosen/software/ClassifyCNV/ClassifyCNV.py'
sys.path.append("/ifs1/home/luohaosen/software/ClassifyCNV/")


def cnv_anno(bed_file, ref_genome_version, out_file):
    # Path(out_file).parent.mkdir(exist_ok=True, parents=True)
    # 运行ClassifyCNV.py
    output_dir = Path(out_file).parent.joinpath("classifycnv_anno/")
    output_dir.mkdir(exist_ok=True, parents=True)
    os.system(f"python {classify_cnv_py} --infile {bed_file} --GenomeBuild {ref_genome_version} --outdir {output_dir}")
    # 读取Scoresheet.txt
    output_score_sheet_file = f"{output_dir}/Scoresheet.txt"
    data = pd.read_table(output_score_sheet_file)
    need_columns_list = ['Chromosome', 'Start', 'End', 'Type', 'Classification', 'Total score', 'Known or predicted dosage-sensitive genes']
    result_data = data[need_columns_list].copy()
    result_data.rename(columns={'Known or predicted dosage-sensitive genes': 'dosage_sensitive_genes',
                                'Total score': 'Total_score'}, inplace=True)
    # 导出
    if Path(out_file).suffix == '.xlsx':
        output_xlsx = f"{Path(out_file).parent}/{Path(out_file).stem}.xlsx"
        result_data.to_excel(output_xlsx, index=False)
    else:
        output_tsv = f"{Path(out_file).parent}/{Path(out_file).stem}.tsv"
        result_data.to_csv(output_tsv, index=False, sep='\t')


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--i', dest='bed_file', help="The bed file from cnv_to_bed.py")
    parser.add_option('--g', '--ref_version', dest='ref_genome_version', help="The ref genome version: {'hg19', 'hg38'}")
    parser.add_option('--o', dest='out_file', help="The output cnv anno file")
    options, args = parser.parse_args()
    cnv_anno(bed_file=options.bed_file, ref_genome_version=options.ref_genome_version, out_file=options.out_file)
