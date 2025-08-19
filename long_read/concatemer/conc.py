import pysam
import os

def extract_softclip_lengths(bam_path, output_txt):
    softclip_lengths = []

    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        for read in bamfile:
            if read.cigartuples:
                for op, length in read.cigartuples:
                    # CIGAR op 4 = soft clip
                    if op == 4:
                        softclip_lengths.append(length)

    with open(output_txt, "w") as f:
        for length in softclip_lengths:
            f.write(f"{length}\n")

bam_t120 = "/users/rg/baygun/BlaER_crg/files_long_read/align_bam_files_lyric/ont-Crg-sc_Hv3_0+_BLaER101Rep1.bam"
bam_t0 = "/users/rg/baygun/BlaER_crg/files_long_read/align_bam_files_lyric/ont-Crg-sc_HpreCap_0+_BLaER101Rep1.bam"

# output
out_dir = "/users/rg/baygun/BlaER_crg/long_read/concatemer/result_concatemer"
os.makedirs(out_dir, exist_ok=True)

extract_softclip_lengths(bam_t120, os.path.join(out_dir, "t120.txt"))
extract_softclip_lengths(bam_t0, os.path.join(out_dir, "t0.txt"))
