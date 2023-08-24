import pysam
import sys
PATH_to_bam_file=sys.argv[1]

bam_file = PATH_to_bam_file
def extract_NM_from_bam(bam_file):

    bam = pysam.AlignmentFile(bam_file, "rb")

    # Iterate over all reads and extract MD or NM label information
    NM_list = []
    for read in bam.fetch():
        if read.has_tag("NM"):
            try:
                NM_tag = read.get_tag("NM")
                NM_list.append(NM_tag)
            except KeyError:
                continue

    bam.close()

    return NM_list


NM_list = extract_NM_from_bam(bam_file)
print(NM_list)
