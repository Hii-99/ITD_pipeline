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
    parser.add_argument('-f', '--file_name', required=True, help="the file ID of tumor sample")
    parser.add_argument('-t', '--TCGA_ID', required=True, help="the TCGA ID of sample")
    parser.add_argument('-p', '--pindel_dir', required=True, help="path to the directory of candidate ITD vep annotated vcf result")
    parser.add_argument('-o', '--output_dir', required=False, default="./", help="the output directory, default is the running directory")
    parser.add_argument('-V', '--Version', action='version', version=VERSION, help="print out the Version of this program and exit")
    parser.add_argument('-v', '--verbose', action='store_true', required=False, default=False, help="turn on the verbosity of this program")
    return parser

def find_line(lines: list[str], keyword: str) -> int:
    for i, line in enumerate(lines):
        if keyword in line:
            return i

def getStartEndPos(maf: pd.DataFrame) -> tuple[int]:
    start_pos = []
    end_pos = []
    for _, row in maf.iterrows():
        if row['Variant_Type'] == 'SNV':
            start_pos.append(int(row['Start_Position']))
            end_pos.append(int(row['Start_Position']))
        elif row['Variant_Type'] == 'insertion':
            start_pos.append(int(row['Start_Position']))
            end_pos.append(int(row['Start_Position']) + 1)
        elif row['Variant_Type'] == 'indel':
            start_pos.append(int(row['Start_Position']))
            end_pos.append(int(row['Start_Position']) + 1)
        elif row['Variant_Type'] == 'deletion':
            start_pos.append(int(row['Start_Position']) + 1)
            end_pos.append(int(row['Start_Position']) + len(row['Reference_Allele']))

    return start_pos, end_pos

def newReference_Allele(maf: pd.DataFrame) -> list:
    change_Reference_Allele = []
    for _, row in maf.iterrows():
        if row['Variant_Type'] == 'INS':
            change_Reference_Allele.append('-')
        else:
            change_Reference_Allele.append(row['Reference_Allele'])
    return change_Reference_Allele

def new_Tumor_Seq_Allele2(maf: pd.DataFrame):
    change_Tumor_Seq_Allele2 = []
    for _, row in maf.iterrows():
        if row['Variant_Type'] == 'DEL':
            change_Tumor_Seq_Allele2.append('-')
        else:
            change_Tumor_Seq_Allele2.append(row['Tumor_Seq_Allele2'])
    return change_Tumor_Seq_Allele2


def mostImportantVariant(maf: pd.DataFrame) -> list:
    maf = maf.reset_index(drop=True)
    change_classification = []
    for i in range(len(maf)):
        if len(maf.loc[i, 'Variant_Classification'].split('&')) == 1:
            change_classification.append(maf.loc[i, 'Variant_Classification'].split('&')[0])
        elif len(maf.loc[i, 'Variant_Classification'].split('&')) != 1:
            largest = maf.loc[i, 'Variant_Classification'].split('&')[0]
            for j in range(len(maf.loc[i, 'Variant_Classification'].split('&'))):
                if IMPORT_DICT[maf.loc[i, 'Variant_Classification'].split('&')[j]] < IMPORT_DICT[largest]:
                    largest = maf.loc[i, 'Variant_Classification'].split('&')[j]
            change_classification.append(largest)

    return change_classification

def newVariant_Classification(maf: pd.DataFrame) -> list:
    change_classification = []

    maf = maf.reset_index(drop=True)
    for i in range(len(maf)):
        if maf['Variant_Type'][i] == 'SNP':
            change_classification.append(maf['Variant_Classification'][i])
        elif maf['Variant_Type'][i] == 'DEL':
            if (maf['Variant_Classification'][i] == 'frameshift_variant') or (maf['Variant_Classification'][i] == 'protein_altering_variant'):
                change_classification.append('Frame_Shift_Del')
            elif maf['Variant_Classification'][i] == 'inframe_deletion':
                change_classification.append('In_Frame_Del')
            else:
                change_classification.append(maf['Variant_Classification'][i])
        elif maf['Variant_Type'][i] == 'INS':
            if (maf['Variant_Classification'][i] == 'frameshift_variant') or (maf['Variant_Classification'][i] == 'protein_altering_variant'):
                change_classification.append('Frame_Shift_Ins')
            elif maf['Variant_Classification'][i] == 'inframe_insertion':
                change_classification.append('In_Frame_Ins')
            else:
                change_classification.append(maf['Variant_Classification'][i])
    return change_classification

def newVariant_Classification2(maf: pd.DataFrame) -> list:

    maf = maf.reset_index(drop=True)
    change_classification = []
    for i in range(len(maf)):
        if maf['Variant_Classification'][i] in ['Splice_Site', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Region', 'Translation_Start_Site', 'Missense_Mutation', 'Intron',
                                                     'Silent', 'RNA', "5'UTR", "3'UTR", 'IGR', "5'Flank", "3'Flank", 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins']:
            change_classification.append(maf['Variant_Classification'][i])
        else:
            change_classification.append('Targeted_Region')
    return change_classification



if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    file_name = args.file_name
    TCGA_ID = args.TCGA_ID
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

    data_FORMAT = vcf['INFO'].str.split('|', expand=True) # the second last columns
    data_FORMAT.rename(columns=dict(zip(data_FORMAT.columns, vcf_lines[annotations].split('Format: ')[1].split('">\n')[0].split('|'))), inplace=True)
    data_maf = pd.concat([vcf.drop(columns='INFO'), data_FORMAT], axis=1)

    data_maf.to_csv(path(output_dir, file_name+".vep103.maf"), index=False, sep='\t')

    maftools_maf = data_maf[['SYMBOL', '#CHROM', 'POS', 'REF', 'ALT', 'VARIANT_CLASS', 'Consequence']]
    maftools_maf = maftools_maf.rename(columns=dict(zip(maftools_maf.columns, ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
                                                                   'Variant_Type', 'Variant_Classification'])))
    maftools_maf["Tumor_Sample_Barcode"] = TCGA_ID

    start_pos, end_pos = getStartEndPos(maftools_maf)
    maftools_maf['Start_Change'] = start_pos # new start position
    maftools_maf['End_Position'] = end_pos
    maftools_maf = maftools_maf[['Hugo_Symbol', 'Chromosome', 'Start_Change', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
                          'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode']] # replace start position
    maftools_maf = maftools_maf.rename(columns={'Start_Change':'Start_Position'})  # rename start position

    maftools_maf['Variant_Type'] = maftools_maf['Variant_Type'].replace({'SNV':'SNP', 'insertion':'INS', 'deletion':'DEL', 'indel': 'INS'})

    # change Reference_Allele
    maftools_maf['Reference_Allele'] = newReference_Allele(maftools_maf)

    # change Tumor_Seq_Allele2
    maftools_maf['Tumor_Seq_Allele2'] = new_Tumor_Seq_Allele2(maftools_maf)

    # get most importance Variant_Classification
    maftools_maf['Variant_Classification'] = mostImportantVariant(maftools_maf)

    # rename Variant_Classification
    maftools_maf['Variant_Classification'] = maftools_maf['Variant_Classification'].replace(MAF_VARIANT_MAPPER_DICT)

    # rename Variant_Classification with Variant Type
    maftools_maf['Variant_Classification'] = newVariant_Classification(maftools_maf)
    maftools_maf['Variant_Classification'] = newVariant_Classification2(maftools_maf)

    maftools_maf.to_csv(path(output_dir, file_name + "maftools.maf"), sep='\t', index=False)