# PGS-拷贝数变异注释

基于ClassifyCNV.py的拷贝数变异注释

- 主要脚本在scripts目录下
- 可用于pgs生产流程的cnv注释流程在scripts/pgs_cnv_anno.py，使用前需检查该脚本的CNV_TO_BED_PY和ANNO_CLASSIFY_CNV_PY路径设置是否符合本地情况
- anno_cnv_classify.py脚本的classify_cnv_py路径设置需要根据ClassifyCNV的环境更改
- 由于ClassifyCNV需要下载外部的数据库，可以定期对ClassifyCNV里的Resources内容做更新，更新脚本可用ClassifyCNV的update_clingen.sh和parse_clingen_tsv.py进行手动更新，位置在201服务器:/ifs1/home/luohaosen/software/ClassifyCNV
- ClassifyCNV-官网位置: https://github.com/Genotek/ClassifyCNV
