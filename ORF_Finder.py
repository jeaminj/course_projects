# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 18:07:37 2023
@authors: Cristopher Quintanilla, Hazael Lemus, Jeamin Jung

Practical Computer Concepts for Bioinformatics 
Fall 2023
Final Group Project

ORF Finder: Program outputs .txt file of DNA sequence(s) open reading frames in FASTA format. Input must be DNA sequence file in FASTA format.
"""
import re, os

'''
Author: H. Lemus
--
@Input: Accepts the name of a 'txt' file containing a DNA summary with sequences in fasta format.
Blank spaces are allowed in the sequences, as well as Upper and Lower cases sequences are recognized equally.
@Output: Returns a dictionary with headers as keys, and 'clean' DNA sequences as values.
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
Author: J. Jung
--
Returns the reverse complement of a DNA sequence  
'''  
def reverseComplement(dnaseq):
    reverse_complement_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    # First transcribes the DNA sequence to its complement strand, then reverses and returns the string
    return ''.join([reverse_complement_dict[nuc] for nuc in dnaseq])[::-1]

'''
Author: J. Jung
--
Gets all 6 frame sequences of the DNA sequence(s) extracted from file, returns list of frame sequences 
'''  
def getFramesList(dna_seq):
    frames = []
    sequence_start_tuple = str(dna_seq[:3])
    reverse_start_tuple = (reverseComplement(dna_seq))[:3]

    frame1 = sequence_start_tuple.find(str((sequence_start_tuple)[0]))
    frame2 = frame1 + 1
    frame3 = frame1 + 2 
    frame4 = reverse_start_tuple.find(((reverse_start_tuple)[0]))  
    frame5 = frame4 + 1 
    frame6 = frame4 + 2 
    
    frame_positions = [frame1,frame2,frame3,frame4,frame5,frame6]

    forward_frames = frame_positions[:3]
    reverse_frames = frame_positions[3:]

    #From each frame position, prints the frame sequence separated by codon:
    #Forward frames
    for pos in forward_frames:
        for_codons = []
        for nuc in range(pos, len(dna_seq)-2,3):
            if pos in forward_frames:
                current_codon = dna_seq[nuc:nuc+3]
                for_codons.append(current_codon)
        frames.append(for_codons)
    #Reverse frames
    for pos in reverse_frames:
        rev_codons = []
        for nuc in range(pos, len(dna_seq)-2,3):
            # Reverse frames
            if pos in reverse_frames:
                current_rev_codon = (reverseComplement(dna_seq))[nuc:nuc+3]
                rev_codons.append(current_rev_codon) 
        frames.append(rev_codons)
        
    return frames

'''
Author: H. Lemus
'''    
def divideFrames(listOfframes):
    f1 = listOfframes[0]
    f2 = listOfframes[1]
    f3 = listOfframes[2]
    f4 = listOfframes[3]
    f5 = listOfframes[4]
    f6 = listOfframes[5]
    
    return f1, f2, f3, f4, f5, f6

'''
Author: C. Quintanilla, H. Lemus, J. Jung
--
Getting the user-define ORFs length threshold 
'''    
def getUsrDefOrfLen():
    while True:
        user_input = input("Enter minimum length in bp for ORFs, (any non-integer will set minimum to 50): ")
        try:
            minorflength = int(user_input)
            if minorflength >= 0:
                return minorflength
            else:
                print("Length cannot be negative, defaulted to 50.")
        except ValueError:
            print("Non-integer input detected. Defaulting to 50.")
            return 50
        
        return minorflength

'''
Author: H. Lemus & J. Jung
--
Searches for valid open reading frame(s) within each frame sequence extracted from getFramesList(), returns orfs stored in list
'''
def getORF(frame, min_length):
    orfs = []
    listOrfs = []
    frameLoyalSeq = '0'.join(frame)
    pattern = re.compile(r'ATG0(?:...0).*?(?:TAG|TGA|TAA)')
    orfs = re.findall(pattern, frameLoyalSeq)

    for orf in orfs:
        cleanORF = orf.replace('0', '')
        if len(cleanORF) >= min_length:
            listOrfs.append(cleanORF)

    return listOrfs

'''
Author: H. Lemus
'''    
def getlistOrfFrames(specifSeq, limitlen):
    frame1, frame2, frame3, frame4, frame5, frame6 = divideFrames(getFramesList(specifSeq))
    
    lorf1 = getORF(frame1, limitlen)
    lorf2 = getORF(frame2, limitlen)
    lorf3 = getORF(frame3, limitlen)
    lorf4 = getORF(frame4, limitlen)
    lorf5 = getORF(frame5, limitlen)
    lorf6 = getORF(frame6, limitlen)

    return lorf1, lorf2, lorf3, lorf4, lorf5, lorf6

'''
Author: C. Quintanilla
--
Formats and writes output of header and orfs to file
'''  
def PrintORFS(orf_list, header_info, frame_number, output_file, seqFromDict):
    with open(output_file, 'a') as file:
        
        for i, orf in enumerate(orf_list, start=1):
            if frame_number not in [4, 5, 6]:
                start_pos = seqFromDict.find(orf)
                start_pos = start_pos + 1
            else:
                start_pos = reverseComplement(seqFromDict).find(orf)
                start_pos = -start_pos - 1
            header = f"> {header_info} | Frame = {frame_number} POS = {start_pos} LEN = {len(orf)}"
            file.write(header + '\n')
            
            # Write the ORF codons in groups of three to the file
            codons = [orf[j:j+3] for j in range(0, len(orf), 3)]
            for k in range(0, len(codons), 15):
                file.write(" ".join(codons[k:k+15]) + '\n')

'''
Author: H. Lemus & J. Jung
--
Prompts user to input name of DNA seq file 
'''
def getFileName():
    curr_dir = os.path.dirname(os.path.realpath(__file__))
    fileObtained = False
    print("Please Enter the name of the Fasta File you want to analyze ('file_name.type'): \n")
    while fileObtained == False:
        targetFile = input()
        if targetFile in os.listdir(curr_dir):
            fileObtained = True
            print("\nFile found!\n")
        else:
            print("Please enter a valid file name\n")
    return targetFile

'''
Author: C. Quintanilla & H. Lemus
'''  
def outputFilePer_DnaSeq_InInputFile(mainDictio, limitLength):
    seqHeaderList = mainDictio.keys()
    countOfFiles = 1 
    outputGenericName = "HL_JJ_CQ_output"
    for item in seqHeaderList:
        dictSequence = mainDictio[item]  
        output_fileName = outputGenericName + '_' + item + ".txt"
        countOfFiles += 1
        # Write the ORFs for each frame to the output file using the list of lists frames_list
        frames_list = getlistOrfFrames(dictSequence, limitLength)
        header_info = item  # Get the header from the dictionary
        # This 'empty' printing steps assures that: if there's is a previously existing "HL_JJ_CQ_output.txt" file,
        # It will be emptied before appending the result of the program, to overlapping outputs.
        with open(output_fileName, 'w') as file:
            file.write('')
        for frame_number, orf_list in enumerate(frames_list, start=1):
            PrintORFS(orf_list, header_info, frame_number, output_fileName,dictSequence)

'''
Author: C. Quintanilla & H. Lemus & J. Jung
'''      
def main():
    current_wd = os.getcwd()
    print("Welcome to ORF finder final project from H. Lemus, J. Jung & C. Quintanilla \n")
    print("Fasta files containing multiple sequences are allowed!\n")
    # Prompt the user for the input file name
    input_file_name = getFileName()
        
    myDict = createDictionary(input_file_name)
    userMinOrfLen = getUsrDefOrfLen()
        
        # Output file per Entry in input file   
    outputFilePer_DnaSeq_InInputFile(myDict, userMinOrfLen)
    print(f"The output files should be available in {current_wd}")
main()
