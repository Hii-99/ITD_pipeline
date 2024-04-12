from pyfaidx import Fasta

def get_genome_base(chromosome, position):
    # Load the human genome sequence
    genome = Fasta("path_to_your_human_genome.fa")  # Provide the path to your human genome fasta file

    # Access the sequence at the specified chromosome and position
    sequence = genome[chromosome][position - 1]  # Adjust position to 0-based index

    return sequence

# Example usage
chromosome = "chrX"
position = 12345
base = get_genome_base(chromosome, position)
print("Base at position", position, ":", base)
