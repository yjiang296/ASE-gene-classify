import pysam
import sys
PATH_to_bam_file=sys.argv[1]

bam_file = PATH_to_bam_file
def extract_NM_from_bam(bam_file):
    # 打开BAM文件
    bam = pysam.AlignmentFile(bam_file, "rb")

    # 遍历所有的reads并提取MD或NM标签信息
    NM_list = []
    for read in bam.fetch():
        if read.has_tag("NM"):
            try:
                NM_tag = read.get_tag("NM")
                NM_list.append(NM_tag)
            except KeyError:
                continue

    # 关闭BAM文件
    bam.close()

    return NM_list

# 示例用法

NM_list = extract_NM_from_bam(bam_file)
print(NM_list)
