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

def reverseComplement(dna_seq): # Returns the reverse complement of a DNA sequence
    reverse_complement_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    # First transcribes the DNA sequence to its complement strand, then reverses and returns the string
    return ''.join([reverse_complement_dict[nuc] for nuc in dna_seq])[::-1]

def getFrames(dna_seq):
    frames = []
    start_codon = re.search("ATG", dna_seq)
    reverse_start_codon = re.search("ATG", (reverseComplement(dna_seq)))
    # Returns all 6 frames of a dna sequence into 'frames' array'
    if start_codon and reverse_start_codon:
       frame1 = start_codon.start() # Gets the 'A' int position of the first "ATG" codon in the dna sequence
       frame2 = frame1 + 1 # Gets the 'T' int position of the first "ATG" codon in the dna sequence
       frame3 = frame1 + 2 # Gets the 'G' int position of the first "ATG" codon in the dna sequence
       frame4 = reverse_start_codon.start() # Gets the 'A' int position of the first "ATG" codon in the reverse complement sequence
       frame5 = frame4 + 1 # Gets the 'T' int position of the first "ATG" codon in the reverse complement sequence
       frame6 = frame4 + 2 # Gets the 'G' int position of the first "ATG" codon in the reverse complement sequence
       
       frame_positions = [frame1,frame2,frame3,frame4,frame5,frame6]
       forward_frames = frame_positions[0:3]
       reverse_frames = frame_positions[3:6]

       #From each frame position, prints the frame sequence separated by codon:
       for pos in frame_positions:
            codons = []
            for nuc in range(pos, len(dna_seq)-2,3):
                # Forward frames
                if pos in forward_frames:
                    current_codon = dna_seq[nuc:nuc+3]
                    codons.append(current_codon)
                # Reverse frames
                elif pos in reverse_frames:
                    current_codon = reverseComplement(dna_seq)[nuc:nuc+3]
                    codons.append(current_codon)
                # Quit appending if stop codon in found
                if current_codon in ["TAG", "TAA", "TGA"]:
                    break
            frames.append(codons)
    return frames

def getORFs(frames):
   orfs = []
    # Examines if an actual ORF is present from the frames stored in 'frames' array:
   for sequence in frames:
       print(sequence) # uncomment to see the 6 frames
       print('')
       orf_codons = []
       stop_codons = ['TAG','TGA','TAA']
       ATG_index = None
       # Scans through all codons, searches for first "ATG" and gets value of its index pos if found
       for index, codon in enumerate(sequence):
           if codon == 'ATG':
               ATG_index = index
               break
        # If the sequence of the frame contains a start codon ('ATG):
       if ATG_index is not None:
           stop_index = None
           for index, codon in enumerate(sequence[ATG_index:]):
               if codon in stop_codons:
                   # 
                   stop_index = ATG_index + index
                   break
           # if the frame sequence also contains a stop codon: ORF!
           if stop_index is not None:
               orf_codons.append(sequence[ATG_index:stop_index + 1]) # ORF sequence between START:STOP codon is added to orf_codons array 
           # If frame sequence contains a start codon, but no stop codon: no ORF
           else:
               orf_codons.append("ORF not found in this frame: Start codon is present, but no Stop codon")
       # If frame sequence doesnt contain a start codon: no ORF
       else:
           orf_codons.append("ORF not found in this frame: No start codon present") 
       orfs.append((orf_codons))
   
   # Function returns any valid orfs or message indicating the frame didnt contain any ORFs
   return orfs

'''
End section of functions written by Jeamin Jung
--> Still W.I.P
'''

#frames = []   
#orfs = []
def main():
    
    myDict = createDictionary("sequence.fasta")
    print(myDict)
    # Jeamin ---------- *comment out to not have ORFs printed* 
    key_names = myDict.keys()
    # ORFS FROM FRAMES
    for key in myDict:
        # Getting and printing frames from :
        dna_seq = myDict[key]
        frames = getFrames(dna_seq)
        orf_count = 1
        for entry in frames:
            print(entry)
        # Getting orfs from frames and printing
        if key in key_names:
            print(f"\nAll ORFs from {key}:")
        orfs = getORFs(frames)
        for orf in orfs:
            print(f"[{orf_count}]" , orf)
            orf_count += 1
    
    print("\n----- FRAME 2: -----")
    print(frames[2])
    print("\n ---- ORF from FRAME 2: ----")
    orf2 = (orfs[0])[0]
    print(orf2)
    print("----")
    print(''.join(orf2).replace(',',' '))
    print(len(''.join(orf2).replace(',',' ')))

    # Jeamin ----------
   
   
main()


