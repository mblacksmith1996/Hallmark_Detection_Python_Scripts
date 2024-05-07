import subprocess
import sys
import os
import time
#Dimorphic retrocopies have been identified. I am inquiring if I can generate more specific breakpoints by comparing the filled site and empty site


def complement(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'a':
        return 't'
    if c == 't':
        return 'a'
    if c == 'c':
        return 'g'
    if c == 'g':
        return 'c'
    # If not ACTGactg simply return same character
    return c    

def revcomp(seq):
    c = ''
    seq = seq[::-1] #reverse
    # Note, this good be greatly sped up using list operations
    seq = [complement(i) for i in seq]
    c = ''.join(seq)
    return c
    
def process_line(line):
    """
    Creates a dictionary of the input line
    """
    #print(line)
    if line[4] == "N":
        return ("N/A","N/A")
    
    line_dict = {}
    line_dict["mCanlor_chrom"] = line[0]
    line_dict["mCanlor_start"] = line[1]
    line_dict["mCanlor_end"] = line[2]
    line_dict["sample_name"] = line[3]
    line_dict["sample_chrom"] = line[4]
    line_dict["sample_start"] = int(line[5])
    line_dict["sample_end"] = int(line[6])
    gene_name = line[7]
    line_dict["hit_presence"] = line[8]

    return line_dict, gene_name
    
def pull_record(file_contents):
    """
    Pulls one record from the input file. Records are in 6 line intervals where each line is the same locus but in a seperate dog/wolf.
    """
    locus_dict = {}
    polymorphic = False
    gene_list = []
    lengths = []
    list_of_lines = []
    for i in range(6):
        line = file_contents.readline().split()
        list_of_lines.append(line)
    print(list_of_lines)
    for i in range(len(list_of_lines)):
        current_line = list_of_lines[i]
        if current_line == []:
            if i == 0:
                return "Empty Record. End of File", ""
            else:
                return "Incomplete Record", ""
        dict_of_line, gene_name = process_line(current_line)
        if dict_of_line == "N/A":
            return "N/A", "N/A"
        #if dict_of_line['sample_name'] == "mcanlor":
        #    mCanlor_presence = dict_of_line['hit_presence']
        #dict_of_line['mCanlor_presence'] = mCanlor_presence
        if gene_name !="N/A":
            gene_list.append(gene_name)
        if dict_of_line["hit_presence"] == "YES":
                length = dict_of_line["sample_end"]-dict_of_line["sample_start"]+1
                lengths.append(length)
        locus_dict[dict_of_line["sample_name"]] = dict_of_line
        if dict_of_line["hit_presence"] != "YES" and dict_of_line["hit_presence"] != "DUPLICATE":
            polymorphic = True
    
        #print(current_line)
    #print(locus_dict)
    if len(set(gene_list)) > 1:
        print(gene_list)
        sys.exit("Too many gene name")
    elif "/" in gene_list[0]:
        return "N/A", "N/A"
    else:
        avg_length = sum(lengths)/len(lengths)
        for key in locus_dict.keys():
            locus_dict[key]["gene_name"] = gene_list[0]
            locus_dict[key]["Expected_Length"] = avg_length
                #print(length)#TODO handle a filter for length.
    #print(locus_dict)
    
    return locus_dict, polymorphic

def process_age(age_out,comp):
    age_out = age_out.split("\n")
    for i in range(len(age_out)):
        #print(age_out[i])
        if age_out[i].startswith("Alignment:"):
            #print(age_out[i+1],age_out[i+2])
            #print((age_out[i+1].split()))
            age_out[i+1] = age_out[i+1].replace("]","").replace("[","").replace(","," ")
            age_out[i+2] = age_out[i+2].replace("]","").replace("[","").replace(","," ")
            
            comp["Refined_Filled_start"] = int(age_out[i+1].split()[4])+comp['filled_start']-comp['flank_dist'] #TODO may have an off by 1 error. Check.
            comp["Refined_Filled_end"] = int(age_out[i+1].split()[7])+comp['filled_start']-comp['flank_dist']
            comp["Refined_empty_start"] = int(age_out[i+2].split()[4])+comp['empty_start']-comp['flank_dist']
            comp["Refined_empty_end"] = int(age_out[i+2].split()[7])+comp['empty_start']-comp['flank_dist']
        elif age_out[i].startswith("Identity at breakpoints:"):
            #print(age_out[i+1],age_out[i+2])
            age_out[i+1] = age_out[i+1].replace("]","").replace("[","")
            age_out[i+2] = age_out[i+2].replace("]","").replace("[","")
            if age_out[i+1].split()[3] != "0":
                #handle
                comp['left_TSD_filled'] = [str(int(age_out[i+1].split()[5].split(",")[0])+comp['filled_start']-comp['flank_dist']),str(int(age_out[i+1].split()[5].split(",")[1])+comp['filled_start']-comp['flank_dist'])]
                comp['right_TSD_filled'] = [str(int(age_out[i+1].split()[7].split(",")[0])+comp['filled_start']-comp['flank_dist']),str(int(age_out[i+1].split()[7].split(",")[1])+comp['filled_start']-comp['flank_dist'])]
            else:
                #no TSD
                comp['left_TSD_filled'] = ["N/A","N/A"]
                comp['right_TSD_filled'] = ["N/A","N/A"]
            if age_out[i+2].split()[3] != "0":
                #handle
                comp['left_TSD_empty'] = [int(age_out[i+2].split()[5].split(",")[0])+comp['empty_start']-comp['flank_dist'],int(age_out[i+2].split()[5].split(",")[1])+comp['empty_start']-comp['flank_dist']]
                comp['right_TSD_empty'] = [int(age_out[i+2].split()[7].split(",")[0])+comp['empty_start']-comp['flank_dist'],int(age_out[i+2].split()[7].split(",")[1])+comp['empty_start']-comp['flank_dist']]
            else:
                #no TSD
                comp['left_TSD_empty'] = ["N/A","N/A"]
                comp['right_TSD_empty'] = ["N/A","N/A"]
        else:
            continue
def process_blat(comp, infile):
    with open(infile,'rt') as blat_input:
        match_count = 0
        for line in blat_input.readlines()[5::]:
            match_count = max(int(line.split()[0]),match_count)
            if int(line.split()[0]) == match_count:
                comp["hit_orientation"] = line.split()[8]
    subprocess.run(f"rm {infile}",shell=True)
    return match_count
def run_age(comp):#, gene_dict): #TODO input this as a file
    canines = {"china":"/nfs/turbo/jmkidddata/genomes/China_UNSW_CanFamBas_1.2/ref/China_UNSW_CanFamBas_1.2.fa", \
    "mischka":"/nfs/turbo/jmkidddata/genomes/UU_Cfam_GSD_1.0/ref/UU_Cfam_GSD_1.0.fa", \
    "mcanlor":"/nfs/turbo/jmkidddata/genomes/mCanLor1.2/ref/mCanLor1.2.fa", \
    "nala":"/nfs/turbo/jmkidddata/genomes/Nala_ASM864105v3/ref/Nala_ASM864105v3.fa", \
    "sandy":"/nfs/turbo/jmkidddata/genomes/Sandy_ASM325472v2/ref/Sandy_ASM325472v2.fa", \
    "zoey":"/nfs/turbo/jmkidddata/genomes/zoey/assemblies/2.3/ref/zoey.2.3.fa"}
    
    cDNAs = {"china":"/nfs/turbo/jmkiddscr/anthony-projects/retrocopy_analysis/china_cDNA/china_gene_orientation_all_cDNA.fa", \
    "mischka":"/nfs/turbo/jmkiddscr/anthony-projects/retrocopy_analysis/mischka_cDNA/mischka_gene_orientation_all_cDNA.fa", \
    "mcanlor":"/nfs/turbo/jmkiddscr/anthony-projects/retrocopy_analysis/mCanLor1.2_cDNA/mCanLor1.2_gene_orientation_all_cDNA.fa", \
    "nala":"/nfs/turbo/jmkiddscr/anthony-projects/retrocopy_analysis/nala_cDNA/nala_gene_orientation_all_cDNA.fa", \
    "sandy":"/nfs/turbo/jmkiddscr/anthony-projects/retrocopy_analysis/sandy_cDNA/sandy_gene_orientation_all_cDNA.fa", \
    "zoey":"/nfs/turbo/jmkiddscr/anthony-projects/retrocopy_analysis/zoey_cDNA/zoey_gene_orientation_all_cDNA.fa"}
    comp['flank_dist'] = 5000
    #gtf_genome = "mischka"
    if abs(comp['filled_start']-comp['filled_end']) > 20000 or abs(comp['empty_start']-comp['empty_end']) > 20000:
        comp["hit_orientation"] = "Locus_Too_Long"
        comp["empty_fill"] = "N/A"
        comp["filled_fill"] = "N/A"
        #comp["TSD_seq"] = "N/A"
        comp["Cleavage_Site"] = "N/A"
        comp['left_TSD_filled'] = ["N/A","N/A"]
        comp['right_TSD_filled'] = ["N/A","N/A"]
        comp["Refined_Filled_start"] = "N/A"
        comp["Refined_Filled_end"] = "N/A"
        comp["Refined_empty_start"] = "N/A"
        comp["Refined_empty_end"] = "N/A"
        #print(comp)
        return comp
    
    extract_cmd = f"samtools faidx {canines[comp['filled_sample']]} {comp['filled_chrom']}:{int(comp['filled_start'])+1-comp['flank_dist']}-{int(comp['filled_end'])+comp['flank_dist']} > filled_coords.fa && samtools faidx {canines[comp['empty_sample']]} {comp['empty_chrom']}:{int(comp['empty_start'])+1-comp['flank_dist']}-{int(comp['empty_end'])+comp['flank_dist']} > empty_coords.fa"
    subprocess.run(extract_cmd,shell=True)
    print(extract_cmd)

    age_cmd = f"age_align filled_coords.fa empty_coords.fa -both"
    proc = subprocess.Popen(age_cmd,shell=True,stdout=subprocess.PIPE)
    age_out,err = proc.communicate()
    age_out = age_out.decode()
    print(age_out)
    process_age(age_out,comp)
    #sys.exit()
    print(comp)

    #validate the insertion
    filled_length = comp["Refined_Filled_end"] - comp["Refined_Filled_start"]
    empty_length = abs(comp["Refined_empty_end"] - comp["Refined_empty_start"])
    
    if filled_length > empty_length:
        #length_dif = abs(mCanlor_length-comp["Expected_Length"])
        ins_length = filled_length
        longer = comp["filled_sample"]
        post_extract = f'samtools faidx {canines[longer]} {comp["filled_chrom"]}:{comp["Refined_Filled_start"]}-{comp["Refined_Filled_end"]} > blat_query.txt'
        
    else:
        #length_dif =  abs(sample_length-comp["Expected_Length"])
        ins_length = empty_length
        longer = comp["empty_sample"]
        post_extract = f'samtools faidx {canines[longer]} {comp["empty_chrom"]}:{comp["Refined_empty_start"]}-{comp["Refined_empty_end"]} > blat_query.txt'
    #print(length_dif)
    #print(post_extract)
    gene_extract = f'samtools faidx {cDNAs[longer]} {comp["gene_name"]} > blat_reference.txt'
    subprocess.run(post_extract,shell=True)
    subprocess.run(gene_extract,shell=True)
    subprocess.call(f"blat -minIdentity=85 blat_query.txt blat_reference.txt blat_output.psl",shell=True)

    #sys.exit()
    match_count = process_blat(comp, "blat_output.psl")

                
    print(match_count)
    if match_count > ins_length *.75:
        match = True
    else:
        match = False
    print(match,longer)
    
    if min(filled_length,empty_length) > 20:
        comp["hit_orientation"] = "Complex_Locus"
        comp["empty_fill"] = "N/A"
        comp["filled_fill"] = "N/A"
        #comp["TSD_seq"] = "N/A"
        comp["Cleavage_Site"] = "N/A"
        #print(comp)
        return comp
        
        #sys.exit("Insertion comes with a corresponding deletion of 10bp or larger")
    if max(filled_length,empty_length) < 80 or match == False:
        subprocess.call(f"blat -minIdentity=85 filled_coords.fa blat_reference.txt blat_filled.psl",shell=True)
        match1 = process_blat(comp, "blat_filled.psl")
        
        subprocess.call(f"blat -minIdentity=85 empty_coords.fa blat_reference.txt blat_empty.psl",shell=True)
        match2 = process_blat(comp, "blat_empty.psl")
        print(match1,match2)
        
        if match1 > ins_length*.5 and match2 > ins_length*.5:
            print("Both sites are filled sites")
            comp["Cleavage_Site"] = "N/A" 
            comp["empty_fill"] = True
            comp["filled_fill"] = True
            return comp
        else:
            print("Insertion is 2x smaller than match. This could cause errors in detection.")
            comp["hit_orientation"] = "Complex_Locus"
            comp["empty_fill"] = "N/A"
            comp["filled_fill"] = "N/A"
            #comp["TSD_seq"] = "N/A"
            comp["Cleavage_Site"] = "N/A"
            #print(comp)
            return comp
            
            #sys.exit()
        #sys.exit("No insertion detected")
        
    if match == True and longer == comp['empty_sample']:
        comp["empty_fill"] = True
        comp["filled_fill"] = False
    elif match == True and longer == comp['filled_sample']:
        comp["empty_fill"] = False
        comp["filled_fill"] = True 
    elif match == False: #TODO size cutoff?
        comp["empty_fill"] = False
        comp["filled_fill"] = False 
        
    comp["insertion_genome"] = longer
    #print(comp)
    
    print(comp)
    if comp["empty_fill"] == False and comp["filled_fill"] == False:
        comp["Cleavage_Site"] = "N/A"
        cmd = ""
    elif comp["empty_fill"] == True and comp["filled_fill"] == True:
        comp["Cleavage_Site"] = "N/A"
        cmd = ""
    elif comp['filled_fill'] == True:
        if comp["left_TSD_filled"][0] == "N/A":
            comp["Cleavage_Site"] = "N/A"
            cmd = ""
        else:
            cmd = f"samtools faidx {canines[longer]} {comp['filled_chrom']}:{int(comp['left_TSD_filled'][0])-2}-{int(comp['left_TSD_filled'][0])+4}"
            subprocess.run(cmd,shell=True)
    elif comp['empty_fill'] == True:
        if comp["left_TSD_empty"][0] == "N/A":
            comp["Cleavage_Site"] = "N/A"
            cmd = ""
        else:
            cmd = f"samtools faidx {canines[longer]} {comp['empty_chrom']}:{int(comp['left_TSD_empty'][0])-2}-{int(comp['left_TSD_empty'][0])+4}"
            subprocess.run(cmd,shell=True)
    else:
        sys.exit("Invalid TSD calculation")
        #print(cmd)
    #if cmd != "":
    #    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    #    faidx_out,err = proc.communicate()
    #    comp["TSD_seq"] = faidx_out.decode().split("\n")[1]

        
    if cmd != "":
        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        out,err = proc.communicate()
        out = out.decode()
        out = out.split("\n")[1].rstrip()
        #print(out)
        
        if comp["hit_orientation"] == "+":
            out = revcomp(out)
        #print(out)
        comp["Cleavage_Site"] = out
    else:
        comp["Cleavage_Site"] = "N/A"
    return comp
    #subprocess.run(age_cmd,shell=True)
def create_comparison_dictionary(filled_sample,empty_sample):
    comparison_dict = \
    {"filled_sample":filled_sample["sample_name"],\
    "filled_chrom":filled_sample["sample_chrom"],\
    "filled_start":filled_sample["sample_start"],\
    "filled_end":filled_sample["sample_end"],\
    "empty_sample":empty_sample["sample_name"],\
    "empty_chrom":empty_sample["sample_chrom"],\
    "empty_start":empty_sample["sample_start"],\
    "empty_end":empty_sample["sample_end"],\
    "gene_name":filled_sample["gene_name"],\
    }
    return comparison_dict
    
def process_locus(locus,final_outfile):#, gene_dict):
    print(locus)
    priority = ["mcanlor","mischka","sandy","nala","china","zoey"]
    filled_genome = []
    empty_genome = []
    for key in locus.keys():
        #if key == "mcanlor":
        if locus[key]["hit_presence"] == "YES":
            filled_genome.append(key)
        else:
            empty_genome.append(key)
    #assign filled site
    #for genome in priority:
        #if 
    for genome in priority:
        if genome in empty_genome:
            empty_comparison = genome
            break
    for genome in priority:
        if genome in filled_genome:
            filled_comparison = genome
            break
    print(empty_comparison,filled_comparison)
    print(filled_genome,empty_genome)
    #time.sleep(5)
    #print(comp,'hi')
    #handle first comparison
    permanent_empty = ""
    while permanent_empty == "":
        comp_dict = create_comparison_dictionary(locus[filled_genome[0]],locus[empty_genome[0]])
        print(comp_dict)
        #if comp_dict["gene_name"] == "MARCKSL1":
        #    return ""
        #if comp_dict["gene_name"] == "PSMA1":
        #    return ""
        #if comp_dict["gene_name"] != "MGST3":
        #    return ""
        #if comp_dict["filled_chrom"] != "chr9":
        #    return ""
        comp = run_age(comp_dict)#,gene_dict)
        
        print(comp)
        
        resolved = {}
        #TODO START HERE NEED TO HANDLE THE FILTERING
        if comp["empty_fill"] == True:
            resolved[comp["empty_sample"]] = True
            extract = "\t".join([comp["filled_chrom"],str(comp["Refined_Filled_start"]),str(comp["Refined_Filled_end"]),'False_Negative',comp['gene_name'],comp['filled_sample'],comp['Cleavage_Site']])+"\n"
            print(extract)
            #sys.exit("Filled empty site")
        elif comp['empty_fill'] == False:
            resolved[comp["empty_sample"]] = False
            extract = ""
        else:
            resolved[comp['empty_sample']] = "Unresolved"
            extract = "\t".join([comp["empty_chrom"],str(comp["Refined_empty_start"]),str(comp["Refined_empty_end"]),comp['hit_orientation'],comp['gene_name'],comp['empty_sample'],comp['Cleavage_Site']])+"\n"
        final_outfile.write(extract)
        print(resolved)

        if resolved[empty_genome[0]] == False:
            permanent_empty = empty_genome[0]
        else:# resolved[empty_genome[0]] == True:
            print("No empty locus detected. Returning resolved as is")
            print(resolved)
            #sys.exit()
        #else:
            empty_genome.pop(0)
            if len(empty_genome) == 0:
                return resolved
            continue
        #locus['mcanlor']["mCanlor_chrom"],locus['mcanlor']["mCanlor_start"],locus['mcanlor']["mCanlor_end"],comp["left_TSD_filled"][0],comp['left_TSD_filled'][1],comp["right_TSD_empty_TSD_filled"][0],comp['right_TSD_filled'][1]
        #locus['mcanlor']["mCanlor_chrom"],locus['mcanlor']["mCanlor_start"],locus['mcanlor']["mCanlor_end"],comp["left_TSD_empty"][0],comp['left_TSD_empty'][1],comp["right_TSD_empty"][0],comp['right_TSD_empty'][1]        
        if comp["filled_fill"] == True:
            resolved[comp["filled_sample"]] = True
            extract = "\t".join([comp["filled_chrom"],str(comp["Refined_Filled_start"]),str(comp["Refined_Filled_end"]),comp['hit_orientation'],comp['gene_name'],comp['filled_sample'],comp['Cleavage_Site'],locus['mcanlor']["mCanlor_chrom"],locus['mcanlor']["mCanlor_start"],locus['mcanlor']["mCanlor_end"],comp["left_TSD_filled"][0],comp['left_TSD_filled'][1],comp["right_TSD_filled"][0],comp['right_TSD_filled'][1]])+"\n"
            #final_outfile.write(extract)
            #print(extract)
            #sys.exit()
        elif comp['filled_fill'] == False:
            resolved[comp["filled_sample"]] = False
            sys.exit("Empty filled site")
        else:
            resolved[comp['filled_sample']] = "Unresolved"
            extract = "\t".join([comp["filled_chrom"],str(comp["Refined_Filled_start"]),str(comp["Refined_Filled_end"]),comp['hit_orientation'],comp['gene_name'],comp['filled_sample'],comp['Cleavage_Site'],locus['mcanlor']["mCanlor_chrom"],locus['mcanlor']["mCanlor_start"],locus['mcanlor']["mCanlor_end"],comp["left_TSD_filled"][0],comp['left_TSD_filled'][1],comp["right_TSD_filled"][0],comp['right_TSD_filled'][1]])+"\n"
            #final_outfile.write(extract)
            #sys.exit()
        final_outfile.write(extract)
        
            #resolved[empty_genome[0]] == "Unresolved"
            #get new empty.
            #sys.exit()
    
    if len(filled_genome) == 1:
        return resolved
    else:
        for i in range(len(filled_genome)-1):
            j = i+1
            comp_dict = create_comparison_dictionary(locus[filled_genome[j]],locus[permanent_empty])
            print(comp_dict)
            comp = run_age(comp_dict)
            if comp["filled_fill"] == True:
                resolved[comp["filled_sample"]] = True
                extract = "\t".join([comp["filled_chrom"],str(comp["Refined_Filled_start"]),str(comp["Refined_Filled_end"]),comp['hit_orientation'],comp['gene_name'],comp['filled_sample'],comp['Cleavage_Site'],locus['mcanlor']["mCanlor_chrom"],locus['mcanlor']["mCanlor_start"],locus['mcanlor']["mCanlor_end"],comp["left_TSD_filled"][0],comp['left_TSD_filled'][1],comp["right_TSD_filled"][0],comp['right_TSD_filled'][1]])+"\n"
                final_outfile.write(extract)
            else:
                resolved[comp["filled_sample"]] = False
                print("Complex_Alignment detected")
                extract = "\t".join([comp["filled_chrom"],str(comp["Refined_Filled_start"]),str(comp["Refined_Filled_end"]),comp['hit_orientation'],comp['gene_name'],comp['filled_sample'],comp['Cleavage_Site'],locus['mcanlor']["mCanlor_chrom"],locus['mcanlor']["mCanlor_start"],locus['mcanlor']["mCanlor_end"],comp["left_TSD_filled"][0],comp['left_TSD_filled'][1],comp["right_TSD_filled"][0],comp['right_TSD_filled'][1]])+"\n"
                final_outfile.write(extract)
                #sys.exit("Empty filled site after previous declaration")
                
                
            print(resolved)
    #sys.exit()
    return resolved
    #sys.exit()

input_file = "/nfs/turbo/jmkiddscr/anthony-projects/retrocopy_analysis/blat_results/prefiltered_polymorphics.for_matt.txt"
#input_file = "/nfs/turbo/jmkiddscr/anthony-projects/retrocopy_analysis/blat_results/polymorphic_rcs.for_matt.txt"
infile = open(input_file,'rt')





continue_iterating = True
final_outfile = open("refined_output.txt",'wt')
while continue_iterating:
    comp, poly = pull_record(infile)
    if comp == "Empty Record. End of File":
        continue_iterating = False
        break
    elif comp == "Incomplete Record":
        sys.exit("Incomplete Record")
    elif comp == "N/A":
        continue
    #print(comp, poly)
    if poly == True:
        process_locus(comp,final_outfile)#, gene_dict)
#print(comp, poly)
sys.exit()

# coordinate_genome = "mcanlor"
# for locus in loci:
    # if "6143660" not in locus:
        # continue
    # locus[4] = locus[4].split(",")
    # print(locus[4])
    # empty_site_canine = []
    # filled_site_canine = []
    # for dog in canines.keys():
        # if dog in locus[4]:
            # filled_site_canine.append(dog)
        # else:
            # empty_site_canine.append(dog)
    # print(empty_site_canine,filled_site_canine)
    
    # #extract flank from coordinate_genome
    # with open("test.fa",'wt') as outfile:
        # cmd = f"samtools faidx {canines[coordinate_genome][0]} {locus[0]}:{locus[2]}-{int(locus[2])+2000}"
        # print(cmd)
        # proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        # faidx,err = proc.communicate()
        # faidx = faidx.decode()
        # for content in faidx.split("\n"):
            # if content.startswith(">"):
                # content = content + f"_{locus[0]}_{locus[2]}_{int(locus[2])+2000}_{coordinate_genome}\n"
            # print(content)
            # outfile.write(content)
    # ooc_file = "/".join(canines[filled_site_canine[0]][0].split("/")[0:-1]) + "/11.ooc"
    # blat_cmd = f"blat {canines[filled_site_canine[0]][0]} test.fa -ooc={ooc_file} filled_test.psl"
    # print(blat_cmd)
    # subprocess.run(blat_cmd,shell=True)
    
    # ooc_file = "/".join(canines[empty_site_canine[0]][0].split("/")[0:-1]) + "/11.ooc"
    # blat_cmd = f"blat {canines[empty_site_canine[0]][0]} test.fa -ooc={ooc_file} empty_test.psl"
    # print(blat_cmd)
    # subprocess.run(blat_cmd,shell=True)
    
    
    # #process the files:
    # total_hits = 0
    # with open("filled_test.psl",'rt') as infile:
        # for line in infile.readlines()[5::]:
            # line = line.split()
            # if line[13] != locus[0] or int(line[0]) < 800:
                # continue
            # total_hits +=1
            # actual_hit_filled = line
    # if total_hits > 1:
        # sys.exit("Too many filled hits")
        
    # total_hits = 0
    # with open("empty_test.psl",'rt') as infile:
        # for line in infile.readlines()[5::]:
            # line = line.split()
            # if line[13] != locus[0] or int(line[0]) < 800:
                # continue
            # total_hits +=1
            # actual_hit_empty = line
    # if total_hits > 1:
        # sys.exit("Too many empty hits")
        
    # print(actual_hit_empty,actual_hit_filled)
    
    # #parse the coordinates
    # hits = [actual_hit_empty,actual_hit_filled]
    # dist = 15000
    # for i in range(len(hits)):
        # hit = hits[i].
        # if hit[8] == "+":
            # coordinates = f"{hit[13]}:{int(hit[15])-dist}-{int(hit[15])+dist}"
        # else:
            # coordinates = f"{hit[13]}:{int(hit[16])-dist}-{int(hit[16])+dist}"
        # if i == 0:
            # empty_fasta_call = f"samtools faidx {canines[empty_site_canine[0]][0]} {coordinates} > Empty_test.fa"
        # else:
            # filled_fasta_call = f"samtools faidx {canines[filled_site_canine[0]][0]} {coordinates} > Filled_test.fa"
    # print(empty_fasta_call,filled_fasta_call)
    # subprocess.run(empty_fasta_call,shell=True)
    # subprocess.run(filled_fasta_call,shell=True)
    
    # if actual_hit_filled[8] != actual_hit_empty[8]:
        # reverse = "-revcom2"
    # else:
        # reverse = ""
    # age_cmd = f"age_align Filled_test.fa Empty_test.fa {reverse}"
    # subprocess.run(age_cmd,shell=True)
    # sys.exit()