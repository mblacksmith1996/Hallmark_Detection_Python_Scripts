import sys
import os
import subprocess
generic_scripts_path = os.path.dirname(__file__)
sys.path.append(generic_scripts_path)
import Detecting_Hallmarks_Functions as Hallmarks
from optparse import OptionParser 


if __name__ == "__main__":
    USAGE = """ USAGE:  python Find_Hallmarks.py --locus_file <str path to input file> --locus <str locus of interest chr:start-end>\
--bam <bam file> --reference <refence fasta file>"""
    
    parser = OptionParser(USAGE)
    parser.add_option("--locus_file",dest='retrogene_file',help='path to file containing loci of interest. 1-based closed is expected')
    parser.add_option("--reference",dest='ref',help='reference fasta file used for minimap2 alignment step')
    (options, args)=parser.parse_args()
    
    proc = subprocess.Popen("module list",shell=True,stdout=subprocess.PIPE)
    proc,err = proc.communicate()
    proc = proc.decode()
    print(proc)   

    #create dict of chrom lengths
    chr_dir = {}
    with open(f"{options.ref}.fai",'rt') as infile:
        for line in infile:
            line = line.split()
            chr_dir[line[0]] = int(line[1])

    #define variables
    dist_from_start = 100
    dist_list = [10,20,30,40,50,60,70,80,90,100,200,300,400,500,750,1000,2000]

    for i in range(len(dist_list)):
        end_of_chrs = 0
        valid_retro = 0
        dist = dist_list[i]
        print(f"Dist={dist}")
        with open(options.retrogene_file,'rt') as input_file:
            with open(f"Retrogene_Output_dist_{dist}.txt",'wt') as outfile:
                #split line and determine orientation of retrogene
                for line in input_file.readlines():
                    print(line.rstrip())
                    split_line = line.split()
                    if split_line[3] == "+":
                        orientation = "Forward"
                    elif split_line[3] == "-":
                        orientation = "Reverse"
                    else:
                        continue
                    if split_line[0] not in chr_dir.keys():
                        sys.exit("Chromosome of retrogene not found in fasta index")
                    else:
                        chr_length = chr_dir[split_line[0]]

                    if int(split_line[1]) <= dist_from_start or chr_length - int(split_line[2]) < dist_from_start:
                        end_of_chrs +=1
                        continue
                    else:
                        valid_retro+=1
                        
                    extract_seq = [split_line[0],int(split_line[1]),int(split_line[2])]    
                    seq1,seq2 = Hallmarks.run_water(dist, dist, extract_seq, options.ref, options.ref,orientation,False,chr_length)
                    print(f"seqs are {seq1} {seq2}")
                    if seq1 != "" and len(seq1[2].replace("-","")) >=5 and len(seq2[2].replace("-","")) >=5:
                        TSDs = [f'{extract_seq[0]}:{int(seq1[0].split("-")[0])+int(seq1[1])-1}-{int(seq1[0].split("-")[0])+int(seq1[3])-1}',seq1[2].replace("-",""),f'{extract_seq[0]}:{int(seq2[0].split("-")[0])+int(seq2[1])-1}-{int(seq2[0].split("-")[0])+int(seq2[3])-1}',seq2[2].replace("-","")]
                        print(TSDs)
                    else:
                        TSDs = ["N/A","N/A","N/A","N/A"]
                        print(TSDs)
                    #detect Poly(A). Poly_A cannot be detected outside of TSD if one exists.
                    if "N/A" in TSDs:
                        poly_a, terminated = Hallmarks.Detect_Poly_As(["N/A"],options.ref,extract_seq,dist,dist,orientation,False,chr_length)
                    else:
                        if orientation == "Forward":
                            poly_a, terminated = Hallmarks.Detect_Poly_As(["N/A"],options.ref,extract_seq,dist,dist,orientation,False,chr_length,int(TSDs[2].split(":")[1].split("-")[0])-1)
                        elif orientation == "Reverse":
                            poly_a, terminated = Hallmarks.Detect_Poly_As(["N/A"],options.ref,extract_seq,dist,dist,orientation,False,chr_length,int(TSDs[0].split(":")[1].split("-")[1])+1)
                    print(f"Printing Detected PolyAs {poly_a}")
                    
                    #Filter poly_A
                    final_poly = []
                    if len(poly_a) >= 1:
                        if orientation == "Forward":
                            if "N/A" in TSDs:
                                final_poly = poly_a[0]
                            else:
                                final_poly = poly_a[-1]
                        else:
                            if "N/A" in TSDs:
                                final_poly = poly_a[-1]
                            else:
                                final_poly = poly_a[0]
                        print(final_poly)
                        poly_coords = f"{split_line[0]}:{final_poly[0]}-{final_poly[2]}"
                        poly_seq = final_poly[1]

                    else:
                        poly_coords = "N/A"
                        poly_seq = "N/A"
                        final_poly = "N/A"
                    if final_poly == []:
                        sys.exit()
                    print(final_poly)
                    
                    if "N/A" not in TSDs:
                        endo_site = Hallmarks.identify_endo_site(orientation,TSDs,options.ref)
                        #sys.exit()
                    else:
                        endo_site = "N/A"
                    
                    #Report Results
                    split_line.append(poly_coords)
                    split_line.append(poly_seq)
                    for item in TSDs:
                        split_line.append(str(item))
                    split_line.append(str(terminated))
                    split_line.append(endo_site)
                    output_line = "\t".join(split_line) + "\n"
                    print(output_line,flush=True)
                    outfile.write(output_line)
                    outfile.flush()
        print(f"Number of loci within {dist_from_start} of the start or end of the chromosome is {end_of_chrs}.")
        