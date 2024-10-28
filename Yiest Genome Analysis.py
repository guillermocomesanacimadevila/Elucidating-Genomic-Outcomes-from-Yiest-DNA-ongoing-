import os
def read_file(file_name):
    with open(os.path.expanduser(file_name), 'r') as f:
        return f.read()


codon_aa_dictionary = {"AUG": "M", "UUU": "F", "UUC": "F", 'UUA': 'L', 'UUG': 'L', 'UCU': 'S',
                                'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
                                'CUC': 'L', 'AUC': 'I', 'GUC': 'V', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
                                'CUG': 'L', 'GUG': 'V', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A', 'CCC': 'P',
                                'ACC': 'T', 'GCC': 'A', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'CCG': 'P',
                                'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
                                'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 'UAA': '*', 'CAA': 'Q',
                                'AAA': 'K', 'GAA': 'E', 'UAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                                'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R',
                                'AGC': 'S', 'GGC': 'G', 'UGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                                'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
                                }

def transcription(seq):
    return seq.replace("T", "U")


def translation(seq):
    amino_acid_sequence = ""

    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon in codon_aa_dictionary:
            amino_acid_sequence += codon_aa_dictionary[codon]
        else:
            amino_acid_sequence += codon

    return amino_acid_sequence


FASTA_file_open = read_file("~/Desktop/BB50242_2022_APC1.2_I_cds.fasta")
RNA_FASTA_file = transcription(FASTA_file_open)

##########################
# Set of required built in functions for genome analysis
# Open Fasta file from "read_file" function
# Could set the dictionary for RNA to AA outside both functions


def counting_genes(seq): # Outputs 20K genes, but human = 20K - so needs fixing
    gene_counter = 0
    i = 0
    stop_codons = ["UAA", "UAG", "UGA"]

    while i < len(seq) - 2:
        codon = seq[i:i + 3]
        if codon == "AUG":
            for j in range(i + 3, len(seq) - 2, 3):
                stop_codon = seq[j:j + 3]
                if stop_codon in stop_codons:
                    gene_counter += 1
                    i = j + 3
        else:
            i += 3

    return gene_counter


def protein_counter(seq):
    protein_counter = 0

    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon in codon_aa_dictionary:
            protein_counter += 1
        else:
            i += 3

    return protein_counter


#def null_gene_counter(seq):
 #   null_gene_counter = 0
  #  stop_codons = ["UAA", "UAG", "UGA"]

   # for i in range(0, len(seq), 3):
    #    codon = seq[i:i + 3]
     #   if codon != "AUG" and not codon.endswith(stop_codons):
      #      null_gene_counter += 1

    #return null_gene_counter


#Counting Leucine encoding codons (CUG vs rest)
print([key for key, value in codon_aa_dictionary.items() if value == "L" and key != "CUG"])
# Prints out ['UUA', 'UUG', 'CUU', 'CUC', 'CUA'] as the remaining leucine encoding codons


#Function to estimate CUG codon frequency within FASTA file
def Leucine_counter(seq): # FIX THIS FUNCTION SO THAT IT OUTPUTS THE REAL COUNT OF CUGs and REMAINING LEUCINE ENDOCING CODONS
    CUG_counter = 0
    Non_CUG_counter = 0
    Leucine_codons_B = ["UUA", "UUG", "CUU", "CUC", "CUA"]


    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon == "CUG":
            CUG_counter += 1
        elif codon in Leucine_codons_B:
            Non_CUG_counter += 1

    return CUG_counter, Non_CUG_counter

AA_encoded_Sequence = translation(RNA_FASTA_file)

def other_leucine_codons(seq):
    Leucine_codons_dict = {"UUA": 0, "UUG": 0, "CUU": 0, "CUC": 0, "CUA": 0}
    CUG_dict = {"CUG": 0}

    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon in Leucine_codons_dict:
            Leucine_codons_dict[codon] += 1
        elif codon in CUG_dict:
            CUG_dict[codon] += 1

    return Leucine_codons_dict, CUG_dict, sum(value for key, value in Leucine_codons_dict.items())


# Need to print gene length (in nucleotides) for all genes with the appropiate starting and ending codons
# print(FASTA_file_open)
print(counting_genes(RNA_FASTA_file))
print(protein_counter(RNA_FASTA_file))
print(Leucine_counter(AA_encoded_Sequence)) # Edit this so that it eliminates all lines that start with a greater than symbol
print(other_leucine_codons(RNA_FASTA_file))

# >lcl|CP034456.1_cds_QBM86630.1_1238 [gene=MPUL0A12750] [locus_tag=METSCH_A12750] [protein=segregation subunit 7] [protein_id=QBM86630.1] [location=3190358..3193078] [gbkey=CDS]