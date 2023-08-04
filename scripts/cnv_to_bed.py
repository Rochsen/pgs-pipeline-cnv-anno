# Time:        2023-07-28
# Author:      HaosenLuo
# FileName:    cnv_to_bed.py
import pandas as pd
import numpy as np
from pathlib import Path
import re
from optparse import OptionParser


# 输入格式，接受组合字符串或者分开填写
def __cnv_input(chrom=None, begin_pos=None, end_pos=None, copy_number=None, gender='', cnv_combined_str=None):
    # 输入方式不同，需保证组合组合字符和分开都可以用
    if cnv_combined_str:
        try:
            chrom_tmp, begin_str, end_str, cn_str, gender_value = re.split(r'[-:,]', cnv_combined_str)
        except ValueError:
            raise Exception(f"{cnv_combined_str} is not full pooling")
    else:
        chrom_tmp, begin_str, end_str, cn_str, gender_value =\
            str(chrom), str(begin_pos), str(end_pos), str(copy_number), gender
    # 染色体质控
    if 'chr' in chrom_tmp:
        chrom_str = chrom_tmp
    else:
        chrom_str = 'chr' + chrom_tmp
    return chrom_str, begin_str, end_str, cn_str, gender_value


# 拷贝数转换为DEL和DUP，非阈值范围的返回“-”
def __copy_number_type_standard(chrom='', copy_number='',
                                gender=None,
                                common_chrom_dup_threshold=None, common_chrom_del_threshold=None,
                                man_chrom_dup_threshold=None, man_chrom_del_threshold=None):
    # 常染色体-DUP最小阈值
    dup_min_threshold_common_cnv = float(common_chrom_dup_threshold) if common_chrom_dup_threshold else 2.8
    # 常染色体-DEL最大阈值
    del_max_threshold_common_cnv = float(common_chrom_del_threshold) if common_chrom_del_threshold else 1.2
    # 性染色体(男性的X和Y染色体按单倍体计算，其余的按女性计算)
    if gender:
        if gender == 'women':
            dup_min_threshold_gender_cnv = float(common_chrom_dup_threshold) if common_chrom_dup_threshold else 2.8
            del_max_threshold_gender_cnv = float(common_chrom_del_threshold) if common_chrom_del_threshold else 1.2
        else:
            dup_min_threshold_gender_cnv = float(man_chrom_dup_threshold) if man_chrom_dup_threshold else 1.8
            del_max_threshold_gender_cnv = float(man_chrom_del_threshold) if man_chrom_del_threshold else 0.4
    else:
        raise Exception(f"gender is NoneType, please input gender param")
    # CNV type
    if chrom not in ['chrX', 'chrY']:
        if float(copy_number) in pd.Interval(left=dup_min_threshold_common_cnv, right=np.inf, closed='left'):
            return 'DUP'
        elif float(copy_number) in pd.Interval(left=-np.inf, right=del_max_threshold_common_cnv, closed='right'):
            return 'DEL'
        else:
            return '-'
    else:
        if float(copy_number) in pd.Interval(left=dup_min_threshold_gender_cnv, right=np.inf, closed='left'):
            return 'DUP'
        elif float(copy_number) in pd.Interval(left=-np.inf, right=del_max_threshold_gender_cnv, closed='right'):
            return 'DEL'
        else:
            return '-'


# 根据cnv文件生产
def cnv_to_bed(input_cnv_file, output_bed, gender=None,
               common_chrom_dup_threshold='2.8', common_chrom_del_threshold='1.2',
               man_chrom_dup_threshold='1.8', man_chrom_del_threshold='0.4'):
    # 阈值引用
    common_threshold_dup, common_threshold_del = common_chrom_dup_threshold, common_chrom_del_threshold
    man_threshold_dup, man_threshold_del = man_chrom_dup_threshold, man_chrom_del_threshold
    # 批处理
    chrom_list, begin_list, end_list, cnv_type_list = [], [], [], []
    chrom_list_append, begin_list_append, end_list_append, cnv_type_list_append = chrom_list.append, begin_list.append, end_list.append, cnv_type_list.append
    # cnv文件路径
    with open(input_cnv_file, 'r') as file_input:
        for line_str in file_input.readlines():
            if re.search(re.compile(r"[XY\d+]", re.I), line_str):
                row = line_str.split('\t')
                # 确保前三列是染色体、起始、终止，以及最后一列是拷贝数
                chrom_str, begin_str, end_str, cn_str, gender_value = __cnv_input(chrom=row[0],
                                                                                  begin_pos=row[1],
                                                                                  end_pos=row[2],
                                                                                  copy_number=row[3],
                                                                                  gender=gender)
                chrom_list_append(chrom_str)
                begin_list_append(begin_str)
                end_list_append(end_str)
                cnv_type_list_append(__copy_number_type_standard(chrom=chrom_str, copy_number=cn_str,
                                                                 common_chrom_dup_threshold=common_threshold_dup,
                                                                 common_chrom_del_threshold=common_threshold_del,
                                                                 man_chrom_dup_threshold=man_threshold_dup,
                                                                 man_chrom_del_threshold=man_threshold_del,
                                                                 gender=gender_value))
    zipped_result = zip(chrom_list, begin_list, end_list, cnv_type_list)
    result_data = pd.DataFrame(zipped_result)
    # 导出
    if Path(output_bed).suffix == '.bed':
        Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
        result_data[result_data[result_data.columns[3]] != '-'].to_csv(output_bed, sep='\t', index=False, header=False)
    else:
        raise Exception(f"BED file is necessary for ClassifyCNV.py, please output suffix select .bed")


# 根据输入的cnv内容生成单bed
def input_to_bed(output_bed,
                 chrom=None, begin_pos=None, end_pos=None, copy_number=None, gender=None, cnv_combined_str_list=None,
                 common_chrom_dup_threshold='2.8', common_chrom_del_threshold='1.2',
                 man_chrom_dup_threshold='1.8', man_chrom_del_threshold='0.4'):
    # 单条cnv分染色体、起始、终止、拷贝数、性别输入
    if all((chrom, begin_pos, end_pos, copy_number, gender)):
        chrom_str, begin_str, end_str, cn_str, gender_value = __cnv_input(chrom, begin_pos, end_pos, copy_number,
                                                                          gender)

        cnv_type = __copy_number_type_standard(chrom_str, cn_str, gender_value,
                                               common_chrom_dup_threshold, common_chrom_del_threshold,
                                               man_chrom_dup_threshold, man_chrom_del_threshold)
        if Path(output_bed).suffix == '.bed':
            Path(output_bed).parent.mkdir(exist_ok=True, parents=True)
            with open(output_bed, 'w') as bed_writer:
                if cnv_type != '-':
                    bed_writer.write(f"{chrom_str}\t{begin_str}\t{end_str}\t{cnv_type}")
        else:
            raise Exception(f"BED file is necessary for ClassifyCNV.py, please output suffix select .bed")
    # 单或多组合cnv字符串输入
    # 单组合：chr1,10000,20000,2.112,man
    # 多组合：chr1,10000,20000,2.112,man|chr1,20000,1220000,1.12,man
    elif cnv_combined_str_list and all((chrom, begin_pos, end_pos, copy_number, gender)) is False:
        # 批处理
        chrom_list, begin_list, end_list, cnv_type_list = [], [], [], []
        chrom_list_append, begin_list_append, end_list_append, cnv_type_list_append = chrom_list.append, begin_list.append, end_list.append, cnv_type_list.append
        for cnv_combined_str in cnv_combined_str_list.split('|'):
            chrom_str, begin_str, end_str, cn_str, gender_value = __cnv_input(cnv_combined_str=cnv_combined_str)
            chrom_list_append(chrom_str)
            begin_list_append(begin_str)
            end_list_append(end_str)
            cnv_type_list_append(__copy_number_type_standard(chrom_str, cn_str,
                                                             gender_value,
                                                             common_chrom_dup_threshold, common_chrom_del_threshold,
                                                             man_chrom_dup_threshold, man_chrom_del_threshold))
        zipped_result = zip(chrom_list, begin_list, end_list, cnv_type_list)
        result_data = pd.DataFrame(zipped_result)
        # 导出
        if Path(output_bed).suffix == '.bed':
            Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
            result_data[result_data[result_data.columns[-1]] != '-'].to_csv(output_bed, sep='\t', index=False,
                                                                            header=False)
        else:
            raise Exception(f"BED file is necessary for ClassifyCNV.py, please output suffix select .bed")
    else:
        raise Exception(f"Not input (chrom, begin_pos, end_pos, copy_number, gender) or cnv_combined_str_list fully")


if __name__ == '__main__':
    parser = OptionParser()
    # cnv文件或外部文件的输入方法
    parser.add_option('--i', dest='input_cnv_file', help="Input cnv file to process ")
    # 单条cnv自定义输出所需参数
    parser.add_option('--chr', dest='chrom', default=None, help="The chromosome of CNV region")
    parser.add_option('--begin', dest='begin_pos', default=None, help="The begin site of CNV region")
    parser.add_option('--end', dest='end_pos', default=None, help="The end site of CNV region")
    parser.add_option('--cn', '--copy_number', dest="copy_number", default=None, help="The copy number")
    # 组合单cnv或多cnv格式字符串所需参数
    parser.add_option('--cnv-string', '--cnv_combined_str_list', dest='cnv_combined_str_list', default=None, type=str,
                      help="The combined cnv string")
    # 必要参数
    parser.add_option('--gender', dest='gender', help="Options {man, woman}")
    parser.add_option('--o', dest='output_bed', help="The output bed file absolute path")
    # 可选参数
    parser.add_option('--pc', '--dup-common-threshold', dest='common_chrom_dup_threshold',
                      help="Like dup threshold min-max", default='2.8')
    parser.add_option('--lc', '--del-common-threshold', dest='common_chrom_del_threshold',
                      help="Like del threshold min-max", default='1.2')
    parser.add_option('--pm', '--dup-man-threshold', dest='man_chrom_dup_threshold',
                      help="Like dup threshold min-max for man", default='1.8')
    parser.add_option('--lm', '--del-man-threshold', dest='man_chrom_del_threshold',
                      help="Like del threshold min-max for man", default='0.4')
    options, args = parser.parse_args()
    # 提供两种构建bed文件的方法
    if options.cnv_combined_str_list or all(
            (options.chrom, options.begin_pos, options.end_pos, options.copy_number, options.gender)):
        input_to_bed(output_bed=options.output_bed,
                     chrom=options.chrom,
                     begin_pos=options.begin_pos,
                     end_pos=options.end_pos,
                     copy_number=options.copy_number,
                     gender=options.gender,
                     cnv_combined_str_list=options.cnv_combined_str_list,
                     common_chrom_dup_threshold=options.common_chrom_dup_threshold,
                     common_chrom_del_threshold=options.common_chrom_del_threshold,
                     man_chrom_dup_threshold=options.man_chrom_dup_threshold,
                     man_chrom_del_threshold=options.man_chrom_del_threshold)
    else:
        cnv_to_bed(input_cnv_file=options.input_cnv_file,
                   output_bed=options.output_bed,
                   gender=options.gender,
                   common_chrom_dup_threshold=options.common_chrom_dup_threshold,
                   common_chrom_del_threshold=options.common_chrom_del_threshold,
                   man_chrom_dup_threshold=options.man_chrom_dup_threshold,
                   man_chrom_del_threshold=options.man_chrom_del_threshold)
