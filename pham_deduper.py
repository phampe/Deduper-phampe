#!/usr/bin/env python

### Deduper 10/18/2022

#import bioinfo
import argparse
import itertools
import gzip
import numpy
import matplotlib
import re

def get_args():
    parser=argparse.ArgumentParser(description="This code will parse through a give SAM file and output a SAM file where there are only unique PCR reads")
    parser.add_argument("-f", help= "SAM file that you will be deduping" )
    parser.add_argument("-o", help= "Output Sam File from input" )
    parser.add_argument("-u", help= "File containing list of UMI" )
    #parser.add_argument("-h","--help" help= "read one" )
    return parser.parse_args()

args=get_args() ### gets you your function 

# l = int(args.l)
input_sam_file = (args.f)
output_sam_file = (args.o)
UMI_file = (args.u)

#print(input_file,output_file,UMI_file)
##############################################################################################################
### dictionary for files
### key = chrom_number, start position, UMI , strand
### value =  count

pcr_duplicate_counter_dict = {}

unknown_umi = 0

###############################################################################################################
### getting UMI into a set 

known_umi = set()

with open(UMI_file,"r") as UMI:
    for i in UMI:
        i=i.strip("\n")
        known_umi.add(i)

### going through our given known umi file. we are then adding the umi to a set to refer to. A set is quicker than lists.

#print(known_umi)

##############################################################################################################
### writing out functions 

##function for parsing out which strand the read is from. If it is the minus or plus strand

def parsing_strand(bitwise_flag):
    flag=bitwise_flag
    if ((flag & 16) ==16):
        strand = "minus"
    else:
        strand = "plus"
    return strand

# test = parsing_strand(16)
# print(test)

##Function for CIGAR string

def cigar_string_parsing(cigar_string,strand,start_position):
    if strand == "plus":
        matches = re.findall(r'(\d+)([A-Z]{1})', cigar_string)
        j=0
        for i in matches:

            if j==0:
                if "S" in i:
                    start_position = int(start_position) - int(i[0]) 
                    j+=1
                else:
                    j+=1
            if j>=1:
                break 
                ## this is making sure that we are only accounting for the first S and not the others as we are only worried about the leftmost
        return(start_position)

    elif strand == "minus":
        matches = re.findall(r'(\d+)([A-Z]{1})', cigar_string)
        j=0
        for i in matches:
            if "M" in i:
                start_position = int(start_position)+int(i[0])
                j+=1

            if "S" in i:
                if j!=0:
                    start_position = int(start_position)+int(i[0])
                    j+=1
                else:
                    j+=1

            if "N" in i:
                start_position = int(start_position)+int(i[0])
                j+=1

            if "D" in i:
                start_position = int(start_position)+int(i[0])
                j+=1

            ### we are going through each of the cigar and add that value to our start position We only whant the S
            ### string when its the last one. This code is account for it by adding making sure its not the first one
        return(start_position)

# a = cigar_string_parsing('15S470N56S','plus',5)
# print(a)

##########################################################################################################

with open(output_sam_file,"w") as output_sam:
    with open(input_sam_file, "r") as sam_file:
            for lines in sam_file:

                if lines.startswith("@") :
                    output_sam.write(lines)
                    ### we want to recreate a SAM file so we would need to save adn write the @ onto the file

                else:
                    line = lines.split()
                    qname_parse = line[0].split(":")
                    qname_umi = (qname_parse[len(qname_parse)-1])
                    ### for this we are just grabbing the UMI first. Then checking to see if it is in our known umi list
                    ### if it isn't then we can easily skip it. 
                    ### this part works, cause if it didn't then we would get many unknown umi

                    if qname_umi in known_umi:
                        umi = qname_umi
                        chromosome = line[2]
                        position_start = int(line[3])
                        strand = str(parsing_strand(int(line[1])))
                        adjusted_position = cigar_string_parsing(line[5],strand,position_start)
                        ### we are setting up our dictionary here with our variables. We now have 4 of our position for our dictionary. 
                        ### in this we also have accounted for the cigar string and adjusted our start position 

                        if (chromosome,adjusted_position,umi,strand) in pcr_duplicate_counter_dict:
                            pcr_duplicate_counter_dict[chromosome,adjusted_position,umi,strand]+=1
                            ### if the values are already in the dictonary then we have a duplicate. 
                            
                        else:
                            pcr_duplicate_counter_dict[chromosome,adjusted_position,umi,strand]=1
                            output_sam.write(lines)
                            ### this is a new pcr that we haven't seen yet. Therefore we need to write it out and record it to the dictionary
                            
                    else:
                        unknown_umi+=1

unique_pcr=len(pcr_duplicate_counter_dict)
total_pcr_reads = sum(pcr_duplicate_counter_dict.values())
duplicates = total_pcr_reads - unique_pcr

print(unknown_umi)
print(unique_pcr,total_pcr_reads,duplicates)
