import argparse
import pandas as pd
from natsort import natsort_keygen
from template.banner import wanglab_banner

# Variant_Classification
# vc_nonSyn: 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'Translation_Start_Site', 'Nonsense_Mutation',
# 'Nonstop_Mutation', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation'
VERSION: str="0.1.1"




IMPORT_DICT = {'transcript_ablation':1,
                  'exon_loss_variant':2,
                  'splice_donor_variant':3,
                  'splice_acceptor_variant':4,
                  'stop_gained':5,
                  'frameshift_variant':6,
                  'stop_lost':7,
                  'start_lost':8,
                  'initiator_codon_variant':9,
                  'disruptive_inframe_insertion':10,
                  'disruptive_inframe_deletion':11,
                  'conservative_inframe_insertion':12,
                  'conservative_inframe_deletion':13,
                  'inframe_insertion':14,
                  'inframe_deletion':15,
                  'protein_altering_variant':16,
                  'missense_variant':17,
                  'conservative_missense_variant':18,
                  'rare_amino_acid_variant':19,
                  'transcript_amplification':20,
                  'splice_region_variant':21,
                  'splice_donor_5th_base_variant':22,
                  'splice_donor_region_variant':23,
                  'splice_polypyrimidine_tract_variant':24,
                  'start_retained_variant':25,
                  'stop_retained_variant':26,
                  'synonymous_variant':27,
                  'incomplete_terminal_codon_variant':28,
                  'coding_sequence_variant':29,
                  'mature_miRNA_variant':30,
                  'exon_variant':31,
                  '5_prime_UTR_variant':32,
                  '5_prime_UTR_premature_start_codon_gain_variant':33,
                  '3_prime_UTR_variant':34,
                  'non_coding_exon_variant':35,
                  'non_coding_transcript_exon_variant':36,
                  'non_coding_transcript_variant':37,
                  'nc_transcript_variant':38,
                  'intron_variant':39,
                  'intragenic_variant':40,
                  'INTRAGENIC':41,
                  'NMD_transcript_variant':42,
                  'upstream_gene_variant':43,
                  'downstream_gene_variant':44,
                  'TFBS_ablation':45,
                  'TFBS_amplification':46,
                  'TF_binding_site_variant':47,
                  'regulatory_region_ablation':48,
                  'regulatory_region_amplification':49,
                  'regulatory_region_variant':50,
                  'regulatory_region':51,
                  'feature_elongation':52,
                  'feature_truncation':53,
                  'intergenic_variant':54,
                  'intergenic_region':55}

MAF_VARIANT_MAPPER_DICT = {'splice_acceptor_variant':'Splice_Site',
                                     'splice_donor_variant':'Splice_Site',
                                     'transcript_ablation':'Splice_Site',
                                     'exon_loss_variant':'Splice_Site',
                                     'stop_gained':'Nonsense_Mutation',
                                     'stop_lost':'Nonstop_Mutation',
                                     'splice_region_variant':'Splice_Region',
                                     'initiator_codon_variant':'Translation_Start_Site',
                                     'start_lost':'Translation_Start_Site',
                                     'missense_variant':'Missense_Mutation',
                                     'coding_sequence_variant':'Missense_Mutation',
                                     'conservative_missense_variant':'Missense_Mutation',
                                     'rare_amino_acid_variant':'Missense_Mutation',
                                     'transcript_amplification':'Intron',
                                     'intron_variant':'Intron',
                                     'INTRAGENIC':'Intron',
                                     'intragenic_variant':'Intron',
                                     'incomplete_terminal_codon_variant':'Silent',
                                     'synonymous_variant':'Silent',
                                     'stop_retained_variant':'Silent',
                                     'NMD_transcript_variant':'Silent',
                                     'mature_miRNA_variant':'RNA',
                                     'exon_variant':'RNA',
                                     'non_coding_exon_variant':'RNA',
                                     'non_coding_transcript_exon_variant':'RNA',
                                     'non_coding_transcript_variant':'RNA',
                                     'nc_transcript_variant':'RNA',
                                     '5_prime_UTR_variant':"5'UTR",
                                     '5_prime_UTR_premature_start_codon_gain_variant':"5'UTR",
                                     '3_prime_UTR_variant':"3'UTR",
                                     'TF_binding_site_variant':'IGR',
                                     'regulatory_region_variant':'IGR',
                                     'regulatory_region':'IGR',
                                     'intergenic_variant':'IGR',
                                     'intergenic_region':'IGR',
                                     'upstream_gene_variant':"5'Flank",
                                     'downstream_gene_variant':"3'Flank"}

def path(*args : str) -> str:
    return "/".join(args)

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="vep_vcf2maf.py",description="using pindel SI and TD result to identify ITD candidate")
    parser.add_argument('-p', '--pindel_dir', required=True, help="path to the directory of candidate ITD vep annotated vcf result")
    parser.add_argument('-o', '--output_dir', required=False, default="./", help="the output directory, default is the running directory")
    parser.add_argument('-V', '--Version', action='version', version=VERSION, help="print out the Version of this program and exit")
    parser.add_argument('-v', '--verbose', action='store_true', required=False, default=False, help="turn on the verbosity of this program")
    return parser

def find_line(lines: list[str], string: str) -> str:
    pass

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    file_name = args.file_name
    pindel_dir = args.pindel_dir
    output_dir = args.output_dir
    verbose = args.verbose
 

    if verbose:
        wanglab_banner()
        print(f"""
{'-'*20}
File Name               : {file_name}
Pindel Directory        : {pindel_dir}
Output Directory        : {output_dir}
{'-'*20}\n
""")
        
    vcf_file = path(pindel_dir, file_name + ".pindel.vep.vcf")
    with open(vcf_file) as f:
        vcf_lines = f.readlines()
        f.close()
    
    header = find_line(vcf_lines, "#CHROM")

    annotations = find_line(vcf_lines, "Consequence annotations from Ensembl VEP")
    vcf = pd.read_csv(vcf_file, sep='\t', header=header, dtype=object).sort_values(by=["#CHROM","POS"],key=natsort_keygen(), ascending=True).drop_duplicates(["#CHROM", "POS"], keep="first").reset_index(drop=True)

    
