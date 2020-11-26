#!/usr/bin/env python
import argparse
import itertools
import re
import gzip
import csv
from argparse import RawTextHelpFormatter

# Get Arguments
def get_args():
    parser = argparse.ArgumentParser(description=str("This script takes a SORTED (by chromosome), unzipped SAM file with known UMI's appended to the first column in the format <ColumnInfo>:<UMI>. File must be from single-end reads. Outputs a SAM file without PCR Duplicates." + "\n" + "Example usage (options provide for inputting an unsorted SAM file and writing PCR duplicates to their own file, but not setting the PCR Duplicate option in the bit flag.):" + "\n" + "python ./CHANG_deduper.py -f <SAM_file_path> -u <UMI_index_file_path> -wd 0 -ns 0 -sbo 1"), formatter_class=RawTextHelpFormatter)
    parser.add_argument("-f", "--file", help="SAM file of single-end reads (unzipped)", required=True)
    parser.add_argument("-u", "--umi", help="File containing known UMI's. No header information, each UMI on its own line.", required=True)
    parser.add_argument("-wd", "--write_dups", help="Optional argument to write the found PCR duplicates to a file called '<filename>_PCRDuplicates.txt' in the current working directory. 0 to write duplicates to separate outfile; 1 to ignore. Default is 1.", required=False)
    parser.add_argument("-sbo", "--set_bits_only", help="Rather than remove all PCR duplicates from the file, simply mark them with the 1024 bit in the bit flag column (as per SAM format). PCR Duplicates can still be written out to a separate file using -wd. 0 to output a file still containing PCR Duplicates but with the 1024 bit set in the flag, 1 to ignore. Default is 1.", required=False)
    parser.add_argument("-ns", "--not_sorted", help="Input an unsorted SAM file -- pysam's .sort() function will sort it for you. User must have pysam library installed. SAM file must not be truncated! 0 to input an unsorted SAM file; 1 to ignore. Default is 1.", required=False)
    return parser.parse_args()


write_dups = 1
ns = 1

args = get_args()
input_f = args.file
umi = args.umi
write_dups = args.write_dups
sbo = args.set_bits_only
ns = args.not_sorted

if sbo != None:
    sbo = args.set_bits_only
    sbo = int(sbo)
if write_dups != None:
    write_dups = args.write_dups
    write_dups = int(write_dups)
if ns != None:
    ns = args.not_sorted
    ns = int(ns)


#Define bad file for PCR Duplicates if writing duplicates to file
if write_dups == 0:
    badfile = open("{}_PCRDuplicates.txt".format(input_f), "w")

if ns == 0:
    import pysam
    #Sort input file
    pysam.sort("-o", input_f, "-O", "sam", input_f)

#Initialize empty dictionary
uniquedict = {}

#Input UMI's to dictionary
with open(umi, "r") as umi:
    for line in umi:
        line = line.strip()
        uniquedict[line] = ()

umi.close()

numnondups = 0
numdups = 0
numstrings = 0

true_pos_dict = {}
true_pos = 0
umi_key = ""
deduped = open('{}_deduped'.format(input_f), "w")

chrom_global = ""

with open(input_f, "rt") as sam:
    '''Main Function: Grab all relevant information from the header line (initial position, chromosome, cigar string, UMI, flag). 
    Adjust true position based on directionality (from the flag) and the CIGAR string. PLEASE BE SURE THE FILE IS OF SINGLE-END READS ONLY. UMI's MUST BE KNOWN.''' 
    for line in sam:
        
        numstrings = numstrings + 1
        if not line.startswith("@"):
            
            header = line.split("\t")
            umi_str = header[0].split(":")
            umi_key = umi_str[7]
            chrom = header[1]
            init_pos = int(header[3])
            flag = int(header[4])
            cigar = header[5]
            
            #Pull out cigar string information and put it into a list
            cigarlist =  re.search("^(\d*)([SDNMI]?)(\d*)([S,D,N,M,I]?)(\d*)([S,D,N,M,I]?)(\d*)([S,D,N,M,I]?)(\d*)([S,D,N,M,I]?)(\d*)([S,D,N,M,I]?)(\d*)([S,D,N,M,I]?)(\d*)([S,D,N,M,I]?)(\d*)([S,D,N,M,I]?)$", cigar)
            if cigarlist != None:
                cigarlist = list(cigarlist.groups())
            else:
                break

            
            #Remove unnecessary parts of the cigar string
            cigar = []
            for i in cigarlist:
                if i != '':
                    cigar.append(i)

            #Reset the dictionary if the chromosome changes
            if chrom != chrom_global:
                true_pos_dict = {}
                #print(chrom_global)
            
            
            #Define directionality based on the flag
            if ((flag & 16) == 16):
                direction = 0
                #Direction 0is Reverse

            else:
                direction = 1
                #Direction 1 is Forward



            #Adjust true position based on cigar string and directionality
            if direction == 0: #Direction is Reverse, calculate true position
                
                #print("direction0")
                m = 0
                for i in range(len(cigar)):
                    #print(cigar[i], "= cigar[i] at i =", i)
                    if "M" in cigar[i]:
                        m = m + 1

                    #Identify the index where the first M occurs in the cigar string
                for i in range(len(cigar)):
                    x = 0
                    if cigar[i] == "M":
                        x = i
                        break
                    if x != 0:
                        break
                    
                #Add all numerical values from the elements that are numerical values after (and including) the first M until the end of the cigar string to the 
                #initial position to find true position. Insertions to the reference genome should not be added to the count.
                #
                # Example: 3S68M should yield true position = initial position + 68
                #
                # Example: 3S52M15D32I71M should yield true position = initial position + 52 + 15 + 71. 
                #          WARNING: Insertions to the reference will automatically be included at the first step. They must be taken out in the following step. 
                #          Therefore, the actual count will first = 52 + 15 + 32 + 71 and after the second step will = 52 + 15 + 32 + 71 - 32
                # 
                # Example: 52M16388N19M should yield true position = initial position OF A REVERSE READ + 52 + 16388 + 19 
                    
                count = 0


                    
                #Add up all values starting with (and including) the first M in the variable "count" then add that to initial position to yield true position
                #print((x+1), len(cigar), 1)
                for i in range(x-1,len(cigar),1):
                    if cigar[i].isnumeric() == True:
                        count = count + int(cigar[i])


                
                #If I is in the cigar string, you need to remove it from the count. It has already been added into the count, so you need to subtract that value from the count       
                #Example: 71M15I71M should have a count of 71 + 71, not 71 + 15 + 71
                if "I" in cigar:
                    for i in range(len(cigar)):
                        if cigar[i] == "I":
                            #IndexI = i
                            count = count - (int(cigar[i - 1]))

                true_pos = init_pos + count

            if direction == 1: #Direction is Forward, calculate true position
                #print("direction1")
                
                if cigar[1] == "M": #Example: 71M3S
                    true_pos = init_pos

                elif cigar[1] == "S": #Example: 3S71M
                    true_pos = init_pos - int(cigar[0])
                    #print("\n \n \n AdjustedForSoftClipped \n \n \n")


            #Determine if the same position and UMI have been located previously on the same chromosome indicating a PCR Duplicate:
            if true_pos in true_pos_dict and true_pos_dict[true_pos] == umi_key:
                #print("PCR DUPLICATE", true_pos, umi_key, true_pos_dict[true_pos])
                numdups = numdups + 1
                #If the user specified to set bits only
                if sbo == 0:
                    flag = int(flag) + 1024
                    flag = str(flag)
                    line = ""
                    for i in range(len(header)):
                        if i !=4:
                            line = line + "\t" + header[i]
                        if i == 4:
                            line = line + "\t" + flag

                    deduped.write(line + "\n")
                #If the user specified outputting PCR Duplicates to their own file
                if write_dups == 0:
                    badfile.write(line)
                else:
                    pass
            
            #If the same position and UMI have not been located previously on the same chromosome, add them to the dictionary and write the read out to the _deduped file
            else:
                
                true_pos_dict[true_pos] = umi_key
                deduped.write(line)
                numnondups = numnondups + 1
            #Assign the chromosome for the fully processed line as chrom_global to compare against the next read's chrom
            chrom_global = chrom
#print(numnondups)
#print(numdups)
#print(numstrings)
sam.close()

