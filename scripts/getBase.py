from pyfaidx import Fasta

def get_genome_base(chromosome, position):
    # Load the human genome sequence
    genome = Fasta("/home/hiiluann99/bin/GenomonITDetector38/GRCh38.d1.vd1.fa")  # Provide the path to your human genome fasta file

    # Access the sequence at the specified chromosome and position
    sequence = genome[chromosome][position - 1]  # Adjust position to 0-based index

    return sequence

# Example usage
chr = ["chr1", "chr1", "chr1"]
pos = [2290573, 2303242, 3200983]

for chro, posi in zip(chr,pos):
    base = get_genome_base(chro, posi)
    print("Base at position", posi, ":", base)
