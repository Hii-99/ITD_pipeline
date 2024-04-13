import pandas as pd
from pyfaidx import Fasta
import argparse
from template.banner import wanglab_banner

VERSION: str="0.1.0"
GENOME_PATH = "/home/hiiluann99/bin/GenomonITDetector38/GRCh38.d1.vd1.fa"
GENOME = Fasta(GENOME_PATH)

FAKE_VCF_HEADER = """##fileformat=VCFv4.0
##fileDate=20240230
##source=pindel
##reference=GRCh38.d1.vd1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=PF,Number=1,Type=Integer,Description="The number of samples carry the variant">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=NTLEN,Number=.,Type=Integer,Description="Number of bases inserted in place of deleted code">
##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference depth, how many reads support the reference">
##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth, how many reads support this allele">"""

def path(*args : str) -> str:
    return "/".join(args)

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ITD2.py",description="convert ITD data to vcf")
    parser.add_argument('-g', '--genomonITD_dir', required=True, help="path to the directory of genomon-ITD result")
    parser.add_argument('-o', '--output_dir', required=False, default="./", help="the output directory, default is the running directory")
    parser.add_argument('-V', '--Version', action='version', version=VERSION, help="print out the Version of this program and exit")
    parser.add_argument('-v', '--verbose', action='store_true', required=False, default=False, help="turn on the verbosity of this program")
    return parser

def get_base(chromosome, position):
    return GENOME[chromosome][position - 1]

def get_INFO(row):
    return f"END={row['position2']};SVLEN={row['PDN_length']};SVTYPE=DUP:TANDEM"

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    ITD_dir = args.genomonITD_dir
    output_dir = args.output_dir
    verbose = args.verbose
 
    if verbose:
        wanglab_banner()
        print(f"""
{'-'*20}
Tumor Sample Directory  : {ITD_dir}
Output Directory        : {output_dir}
{'-'*20}\n
""")
    
    ITD_data = pd.read_csv(ITD_dir, index_col = None)
    ITD_data = ITD_data[ITD_data["in_normal"] == 0]

    tumor_type: str = ITD_data["tumor_type"][0]
    samples = ITD_data["sample_id"].unique()

    for sample in samples:
        this_data = ITD_data[ITD_data["sample_id"] == sample]
        this_vcf: pd.DataFrame = this_data[["chromosome1", "position1"]]
        this_vcf.columns = ['#CHROM', 'POS']
        this_vcf["ID"] = '.'
        this_vcf["REF"] = this_vcf.apply(lambda row: get_base(row['#CHROM'], row['POS']), axis=1)
        this_vcf["ALT"] = this_vcf["REF"].astype(str) + this_data["observed_inserted_nucleotide"]*2
        this_vcf["QUAL"] = '.'
        this_vcf["FILTER"] = '.'
        this_vcf["INFO"] = this_data.apply(get_INFO, axis=1)
        this_vcf["FORMAT"] = "."
        with open(path(output_dir,sample + ".filtered.vcf"), 'w') as vcf_file:
            vcf_file.write(FAKE_VCF_HEADER + '\n')
            this_vcf.to_csv(vcf_file, sep='\t', index=False)



