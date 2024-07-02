import os
import subprocess
import sys
import argparse

parser = argparse.ArgumentParser(prog="Pringing_locus.py",description="Prints loci associated with retrocopies. Seperates the TSD from the intervening sequence.")
parser.add_argument("--input_file", dest="locus_file",required=True,help="Path to loci")
parser.add_argument("--reference_path",dest="reference_path",required=True,help="Path to reference containing retrocopies of interest")
args = parser.parse_args()

def split_the_locus(locus):
    locus = locus.replace("-",":").split(":")
    locus[1] = int(locus[1])
    locus[2] = int(locus[2])
    return locus

def launch_samtools_command(reference,locus,detailed):
    cmd = f"samtools faidx {reference} {locus}"
    extract = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    extract,err = extract.communicate()
    extract = extract.decode()
    output = extract.split("\n")
    
    for line in output:
        if line.startswith(">"):
            header = line.rstrip()
            break
    seq = ""
    for line in output:
        if not line.startswith(">"):
            seq = seq + line.rstrip()
    if detailed == True:
        print(header)
    print(seq,'\n')
    #print(output)
  

detailed = True  
with open(args.locus_file,'rt') as infile:
    count = 0
    for line in infile:
        if count > 10:
            break
        print(line)
        line = line.split()
        line[1] = int(line[1])
        line[2] = int(line[2])
        
        #determine orientation
        # if line[6] != "N/A":
            # if line[6].upper() == "AAA":
                # orientation = "FWD"
            # else:
                # orientation = "REV"
        # else:   
            # orientation = "N/A"
            # continue
            
        #create the output
        print(f"Starting Locus {line[0:3]} \n")
        print(line,"\n")
        print("Upstream \n")
        if line[7] != "N/A":
            left_tsd = split_the_locus(line[7])
            right_tsd = split_the_locus(line[9])
            
            launch_samtools_command(args.reference_path,f"{line[0]}:{left_tsd[1]-30}-{left_tsd[1]-1}",detailed)
            print("Left TSD")
            launch_samtools_command(args.reference_path,f"{line[0]}:{left_tsd[1]}-{left_tsd[2]}",detailed)            
            if line[1] - left_tsd[2] > 1:
                print("sequence between the TSD and start of locus")
                launch_samtools_command(args.reference_path,f"{line[0]}:{left_tsd[2]+1}-{line[1]-1}",detailed)            
            print("printing the locus itself")
            launch_samtools_command(args.reference_path,f"{line[0]}:{line[1]}-{line[2]}",detailed)            
            
            if right_tsd[1] - line[2] > 1:
                print("sequence between the end of locus and start of TSD")
                launch_samtools_command(args.reference_path,f"{line[0]}:{line[2]+1}-{right_tsd[1]-1}",detailed)
            print("Right TSD")
            launch_samtools_command(args.reference_path,f"{line[0]}:{right_tsd[1]}-{right_tsd[2]}",detailed)
            
            print("Downstream")
            launch_samtools_command(args.reference_path,f"{line[0]}:{right_tsd[2]+1}-{right_tsd[2]+30}",detailed)
        #print(orientation)
        #sys.exit()
        count +=1