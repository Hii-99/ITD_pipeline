import pandas as pd
import argparse
from template.banner import wanglab_banner

VERSION: str="0.1.0"

VCF_COLNAME: list[str] = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

def path(*args : str) -> str:
    return "/".join(args)

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ITD2.py",description="convert ITD data to vcf")
    parser.add_argument('-g', '--genomonITD_dir', required=True, help="path to the directory of genomon-ITD result")
    parser.add_argument('-o', '--output_dir', required=False, default="./", help="the output directory, default is the running directory")
    parser.add_argument('-V', '--Version', action='version', version=VERSION, help="print out the Version of this program and exit")
    parser.add_argument('-v', '--verbose', action='store_true', required=False, default=False, help="turn on the verbosity of this program")
    return parser

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
Pindel Directory        : {pindel_dir}
Output Directory        : {output_dir}
{'-'*20}\n
""")
    
    ITD_data = pd.read_csv(ITD_dir, index_col = None)
    ITD_data = ITD_data[ITD_data["in_normal" == 0]]

    tumor_type: str = ITD_data["tumor_type"][0]
    samples = ITD_data["sample_id"].unique()

    for sample in samples:
        this_data = ITD_data[ITD_data["sample_id"] == sample]
        this_vcf = this_data[["chromosome1", "position1"]]
        this_vcf.columns = VCF_COLNAME[:2]
        this_vcf["ID"] = 
