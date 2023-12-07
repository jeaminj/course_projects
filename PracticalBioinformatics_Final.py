# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 18:23:08 2023

@author: hazle
"""

import os, re

'''
Reading file Function: createDictionary()
@Input: Accepts the name of a 'txt' file containing a DNA summary with sequences in fasta format.
Blank spaces are allowed in the sequences, as well as Upper and Lower cases sequences are recognized equally.
@Output: Returns a dictionary with headers as keys, and 'clean' DNA sequences as values.

@auth H. Lemus

'''

def createDictionary(fileName):
    # Opening and reading the text file to obtain the totality of its content as a string
    curr_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(curr_dir, fileName), 'r') as myFile:
        content = myFile.read()
    myFile.close()
    
    # The totality of the file content gets split by entries (headers plus DNA)
    # The split function sometimes gives an extra empty space '' at the begining of the file, so we filtered out.
    contEntries = list(filter(None, re.split(r'>', content)))
    keysList = []
    seqList = []
    newDirectory = {}
    
    # Each entry is analyzed to identify its header and DNA sequence separately
    for entry in contEntries:
        
        # Each entry is then split in 2: header and rest of code by detecting the first \n
        linesInEntry = entry.split('\n', 1)
        # The header is orderly stored in a key list
        keysList.append(linesInEntry[0])
        
        # The rest of code is analyzed to find all strings that match a DNA sequence, accounting for blank spaces and case sensitivity
        seqHit = re.findall(r'[ATGCatgc\s]+', linesInEntry[1])
        # All the matches are concatenated to make a single string
        sequence = ''.join(seqHit)
        # The sequence is cleaned
        sequence = re.sub(r'[\s\n]', '', sequence)
        # And an uppercase version of it is stored in a sequence list
        seqList.append(sequence.upper())
            
    # Finally both key / sequence lists are used to create a Directory that is returned.
    newDirectory = dict(zip(keysList, seqList))
    return newDirectory

'''
Example of how to run the function
Should be deleted for final version
'''

'''
Begin section of functions written by Jeamin Jung

'''
def getORF_length():
    # Asks user for minimum ORF length to include, signifies to user the default is 50 bases
    orf_length = input("What is the minimum length of ORFs you wish to return from file? (default = 50 bases): ")
    return orf_length

# TO be utilized for getting reading frames 4,5, and 6
def reverseComplement(dna_seq):
    reverse_complement_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    # First transcribes the DNA sequence to its complement strand, then reverses and returns the string
    return ''.join([reverse_complement_dict[nuc] for nuc in dna_seq])[::-1]

# Orf finder function v1.0 --------------------------------------------------------------------
# *** Currently only prints ONE reading frame, Need to add code for getting other 5 reading frames ***
def findORFs(dna_seq, orf_length):
   start_codon = re.search("ATG", dna_seq)
   # Only begins extrating codons from reading frame if a start codon exists in the dna sequence
   if start_codon:
       # Gets the int position of the first "ATG" codon in the dna sequence
       start_codon_pos = start_codon.start()
       codons = []
       for pos in range(start_codon_pos, len(dna_seq)-2,3):
           current_codon = dna_seq[pos:pos+3]
           codons.append(current_codon)
           # Quit appending if stop codon in found
           if current_codon in ["TAG", "TAA", "TGA"]:
               break 
       # Fulfilling the minimum orf length requirement and retracting output if user-input miminum orf_length exceeds that of actual orf length
       rf = '' 
       joined_rf = rf.join(codons) 
       rf_length = len(joined_rf) 
       if rf_length >= 50:
           if rf_length >= orf_length:
               print(codons)
           else: 
               print(f"The ORF length of this DNA sequence is below your minimum inputted threshold of >{orf_length}< bases")
       else:
           print("The ORF length of this DNA sequence is less than 50 bases")
   else:
       print("No start codon found in this DNA sequence")

'''
End section of functions written by Jeamin Jung
--> Still W.I.P
'''

    
def main():
        
    myDict = createDictionary("TestSequences.txt")

    # Jeamin ----------
    orf_length = int(getORF_length())
    for key in myDict:
        dna_seq = myDict[key]
        findORFs(dna_seq, orf_length)
    # Jeamin ----------
    
main()
