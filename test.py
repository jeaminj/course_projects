import re

def reverseComplement(dna_seq): # Returns the reverse complement of a DNA sequence
    reverse_complement_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    # First transcribes the DNA sequence to its complement strand, then reverses and returns the string
    return ''.join([reverse_complement_dict[nuc] for nuc in dna_seq])[::-1]

dna_seq = 'GGCATGCTT'
sequence_start_tuple = dna_seq[:3]
start_codon = re.search("GGC", dna_seq)
print(sequence_start_tuple)
frame1 = sequence_start_tuple.find(str((sequence_start_tuple)[0]))
frame11 = start_codon.start()
print(frame1)
print(frame11)

reverse_start_tuple = reverseComplement(dna_seq[:3])
print(reverse_start_tuple)
