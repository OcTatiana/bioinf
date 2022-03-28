import sys


standard_starts = ["ATG", "CTG", "TTG"]
standard_stops = ["TAA", "TAG", "TGA"]

mitochondrial_starts = ["ATT", "ATC", "ATA", "ATG", "GTG"]
mitochondrial_stops = ["TAA", "TAG", "AGA", "AGG"]

standard_transl_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

mitochondrial_transl_table = {
    'ATA': 'M', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': '*', 'AGG': '*',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'W', 'TGG': 'W',
}

def find_stop_codon(sequence, stops):
    for i in range(0, len(sequence), 3):
        if sequence[i:i+3] in stops:
            return i+3
    return -1


def orf_finger(sequence, list_of_start_codons, list_of_stop_codons):
    for i in range(len(sequence) - 2):
        if sequence[i:i+3] in list_of_start_codons:
            start_index = i
            stop_index = find_stop_codon(sequence[start_index + 3:], list_of_stop_codons)
            if stop_index > 0:
                return sequence[start_index:stop_index + start_index], stop_index + start_index + 3

    return "don't find", -1


def find_all_orfs(sequence, starts, stops):
    list_of_orfs = []

    i = 0
    while i < len(sequence):
        orf = orf_finger(sequence[i:], starts, stops)

        if orf[0] == "don't find":
            return list_of_orfs
        else:
            list_of_orfs.append(orf[0])
            i += orf[1]

    return list_of_orfs


def nucl_to_aminoacids(nucl_sequence, transl_table):
    prot_sequence = "M"
    for i in range(3, len(nucl_sequence), 3):
        codon = nucl_sequence[i:i + 3]
        prot_sequence += transl_table[codon]
    return prot_sequence


def read_from_fasta(file):
    return file.readlines()[1].replace("\n", "").replace("\r", "")


def write_to_fasta(file, nucl_seq, prot_seq, num):
    file.write(">ORF_" + str(num) + "_nucl\n")
    file.write(nucl_seq + "\n")
    file.write(">ORF_" + str(num) + "_prot\n")
    file.write(prot_seq + "\n")


def main():
    # command line: python main.py name_of_input_fasta_file genome_type
    # Genome_type: standard or mitochondrial

    filename = sys.argv[1]
    fin = open(filename, 'r')
    fout = open('output.fa', 'w')

    sequence = read_from_fasta(fin)
    genome_type = sys.argv[2]

    if genome_type == "standard":
        list_of_orf = find_all_orfs(sequence, standard_starts, standard_stops)
        num = 1
        for orf in list_of_orf:
            write_to_fasta(fout, orf, nucl_to_aminoacids(orf, standard_transl_table), num)
            num += 1
    elif genome_type == "mitochondrial":
        list_of_orf = find_all_orfs(sequence, mitochondrial_starts, mitochondrial_stops)
        for orf in list_of_orf:
            write_to_fasta(fout, orf, nucl_to_aminoacids(orf, mitochondrial_transl_table), 1)
    else:
        print("There isn't such a genome type")


if __name__ == '__main__':
    main()
