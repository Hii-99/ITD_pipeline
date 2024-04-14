import pandas as pd
import numpy as np
import os
import io
import re
import argparse
from template.banner import wanglab_banner

VERSION: str="0.1.1"
THRESHOLD: int = 5

def path(*args : str) -> str:
    return "/".join(args)

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="normal_filtering.py",description="using pooled normal candidate to filter tumor specific ITD")
    parser.add_argument('-f', '--file_name', required=True, help="the file ID of tumor sample")
    parser.add_argument('-g', '--genomonITD_dir', required=True, help="path to the directory of genomon-ITD result")
    parser.add_argument('-n', '--normal_candidate', required=True, help="path to the Normal Candidate result")
    parser.add_argument('-o', '--output_dir', required=False, default="./", help="the output directory, default is the running directory")
    parser.add_argument('-V', '--Version', action='version', version=VERSION, help="print out the Version of this program and exit")
    parser.add_argument('-v', '--verbose', action='store_true', required=False, default=False, help="turn on the verbosity of this program")
    return parser

def split_normal(Normal_candidate: pd.DataFrame) -> dict[pd.DataFrame]:
    NormalDict = dict()
    for chromosome in set(Normal_candidate.chromosome1):
           this_data = Normal_candidate[Normal_candidate.chromosome1 == chromosome]
           NormalDict[chromosome]=this_data
    return NormalDict

def check_normal(row, NormalDict: dict[pd.DataFrame]) ->int:
    chromo = row["chromosome1"]
    start = row["position1"]
    end = row["position2"]
    in_normal = 0

    if chromo in NormalDict:
        thisChrNormal:pd.DataFrame = NormalDict[chromo]
    else:
        return 0

    for _, row in thisChrNormal.iterrows():
        if (abs(start - row["position1"]) <= THRESHOLD) and (abs(end - row["position2"]) <= THRESHOLD):
            in_normal = 1
            break

    return in_normal

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    file_name = args.file_name
    tumor_dir = args.genomonITD_dir
    normal_candidate = args.normal_candidate
    output_dir = args.output_dir
    verbose = args.verbose
 

    if verbose:
        wanglab_banner()
        print(f"""
{'-'*20}
File Name               : {file_name}
Tumor Sample Directory  : {tumor_dir}
Normal Candidate        : {normal_candidate}
Output Directory        : {output_dir}
{'-'*20}\n
""")
        
    g_data = pd.read_csv(path(tumor_dir,file_name+"_tidy.csv"), header=0, index_col=None)
    Normal_dict = split_normal(pd.read_csv(normal_candidate, index_col=None))
    g_data["in_normal"] = g_data.apply(check_normal, axis=1, args=(Normal_dict,))
    g_data.to_csv(path(output_dir, file_name+".filtered.csv"))
