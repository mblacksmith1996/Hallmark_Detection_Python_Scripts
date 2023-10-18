import sys
import os
import subprocess

def process_file(input_file):
    """
    Extract the query, reference, and between sequences from a smith-waterman alignment
    
    Args:
        input_file (str): Path to file from which the alignment is being extracted
    
    Returns:
        seq1: (list) containing the first 13 characters of the a sequence, start position the aligned portion of the A sequence, the aligned nucleotides of the A sequence, the end position of the aligned A sequence
        seq2: (list) containing the first 13 characters of the b sequence, start position the aligned portion of the B sequence, the aligned nucleotides of the B sequence, the end position of the aligned B sequence
    """
    with open(input_file) as infile:
                file_contents = infile.readlines()[32:-3:]
                i = 0
                seq1 = ["","","",""]
                seq2 = ["","","",""]
                while i < len(file_contents):
                    content = file_contents[i].split()
                    if i%4 == 0:
                        if seq1[0] == "":
                            seq1[0] = (content[0])
                            seq1[1] = (content[1])
                        seq1[2] = seq1[2] + content[2]
                        if i+4 == len(file_contents):
                            seq1[3] = (content[3])
                    elif i%4 == 2:
                        if seq2[0] == "":
                            seq2[0] = (content[0])
                            seq2[1] = (content[1])
                        seq2[2] = seq2[2] + content[2]
                        if i+2 == len(file_contents):
                            seq2[3] = (content[3])
                    i+=1
    return(seq1,seq2)

def investigate_TSD_validity(flank_dist,dist,seq1,seq2,length,search_internal,interior_dist):
    """
    Detects if target site duplications are near within a predetermined distance of the locus of interest.
    Additionally, determines if flanks need to be extended (this occurs when the TSD is within the predetermined distance, but can be extended beyond the currently extracted sequence
    
    Args:
        flank_dist: (int) the number of base pairs that are extracted beyond the locus of interest 
        dist: (int) detected TSDs must be within dist bp of the start of the locus to be considered valid_TSD
        seq1: (list) containing the first 13 characters of the a sequence, start position the aligned portion of the A sequence, the aligned nucleotides of the A sequence, the end position of the aligned A sequence
        seq2: (list) containing the first 13 characters of the a sequence, start position the aligned portion of the B sequence, the aligned nucleotides of the B sequence, the end position of the aligned B sequence
        length: (int) length of the insertion. This will be used in future version in which non-reference hallmark detection is possible
        search_internal: (bool) Determines if TSDs can be detected within the locus of interest. Should be set to False for retrogene detection. Later, True will be added for non-reference hallmark detection
        interior_dist: (int) The number of base pairs within the locus of interest that will be used to search for TSDs. Can be used when search_internal is False
        
    Returns:
        bool: True if the flank has been extended and the alignment is to be re-performed.
        bool: True if one or both of the sequences are invalid (as determined by proximity to the locus of interest). If True, the above bool is also set to False     
    """
    alignment_fail = []

    if search_internal == True:
        for seq in [seq1,seq2]:
            if abs(flank_dist - int(seq[1])) > dist and abs(flank_dist - int(seq[3])) > dist:
                if seq == seq1:
                    alignment_fail.append("seq1")
                else:
                    alignment_fail.append("seq2")
    else:
        print(seq1,seq2)
        print(dist,flank_dist)
        if abs(flank_dist - int(seq1[3])) > dist and abs(flank_dist - int(seq1[1])) > dist:
            alignment_fail.append(seq2)
        if abs(interior_dist - int(seq2[1])) > dist and abs(interior_dist - int(seq2[3])) > dist:   
            alignment_fail.append(seq1)
    print(alignment_fail)

    if len(alignment_fail) >= 1:
        print(f"At least 1 alignment failed to be withing {dist}bp of the start/end of the inserted sequence")
        print(seq1,seq2)
        return False, False
        
    else:
        if search_internal == True:
            if (int(seq1[1]) < 5) or (int(seq2[1]) < 5) or (flank_dist*2 - int(seq1[3]) < 5) or (flank_dist*2 - int(seq2[3]) < 5):
                #print(length)
                if flank_dist > length/2:
                    False,False
                    #sys.exit()
                #print(seq1,seq2)
                print("TSD may be longer than extracted flank. Doubling flank and trying again")
                subprocess.run("rm water.txt",shell=True)
                return True, False
        else:
            #sys.exit()
            if (int(seq1[1]) < 5) or (flank_dist+interior_dist - int(seq2[3]) < 5):
                print("TSD may be longer than extracted flank. Extending flank and trying again")
                subprocess.run("rm water.txt",shell=True)
                return True, False
        #return False, True
        
        
def run_water(flank_dist, dist, coords, extraction_path,ref,orientation, search_internal, chrom_length):
    """
    Function to extract two sequences from a fasta file of interest and then perform a smith-waterman alignment using the EMBOSS module
    
    Args:
        flank_dist: (int) the number of base pairs that are extracted beyond the locus of interest 
        dist: (int) detected TSDs must be within dist bp of the start of the locus to be considered valid_TSD
        coords: (list) containing the sequence name, start base, and end base of the locus of interest. 1 based 
        extraction_path: (str) path to the fasta containing the sequence being extracted
        ref: (str) path to reference. As of now not used, but likely will be for non-reference detection
        orientation: (str) Forward or Reverse, the orientation of the retroelement in reference to the chromosome or contig
        search_internal: (bool) Determines if TSDs can be detected within the locus of interest. Should be set to False for retrogene detection. Later, True will be added for retrogene hallmark detection
        chrom_length: the length of the chromosome or contig in which the locus of interest resides
        
    Returns:
        seq1: (list) containing the first 13 characters of the a sequence, start position the aligned portion of the A sequence, the aligned nucleotides of the A sequence, the end position of the aligned A sequence, returns "" if no TSD detected
        seq2: (list) containing the first 13 characters of the b sequence, start position the aligned portion of the B sequence, the aligned nucleotides of the B sequence, the end position of the aligned B sequence, returns "" if no TSD detected
    """
    interior_dist = 0
    re_run = True
    while re_run:
        if search_internal == True:
            extract_up = f"{coords[0]}:{max(coords[1]-flank_dist,1)}-{coords[1]+(flank_dist -1)}"
            extract_down = f"{coords[0]}:{coords[2]-(flank_dist-1)}-{coords[2]+flank_dist}"
        elif search_internal != True:
            dist = flank_dist
            extract_up = f"{coords[0]}:{max(coords[1]-dist,1)}-{coords[1]+(interior_dist-1)}"
            extract_down = f"{coords[0]}:{coords[2]-(interior_dist-1)}-{coords[2]+dist}"

        command = f"samtools faidx {extraction_path} {extract_up} > upstream_and_5_start.fa ; samtools faidx {extraction_path} {extract_down} > downstream_and_3_end.fa"
        print(command)
        subprocess.run(command, shell=True)    
        gap_open = 10
        gap_extend = 10
        custom_matrix_path = os.path.dirname(__file__)
        cmd = f"water upstream_and_5_start.fa downstream_and_3_end.fa -gapopen {gap_open} --gapextend {gap_extend} -datafile {custom_matrix_path}/Custom_Matrix -outfile water.txt"
        #print(cmd)
        subprocess.run(cmd, shell=True)
        seq1,seq2 = process_file("water.txt")   
        re_run, valid_TSD = investigate_TSD_validity(flank_dist,dist,seq1,seq2,coords[2]-coords[1]+1,search_internal,interior_dist)
        
        #Check for edge cases
        if coords[1]-flank_dist <= 1:
            re_run = False
        elif coords[2]+dist >= chrom_length:
            re_run = False
        elif len(seq1[2]) == 1:
            re_run = False
        elif "N" in seq1[2] or "N" in seq2[2]:
            sys.exit("Masked based in aligned sequence.")
        if re_run == True:
            flank_dist = flank_dist+5
            continue
            
        if valid_TSD != True:
            print("No Valid TSD detected")
            return "",""
        else:
            #print(seq1,seq2)
            if seq1[2].lower() == seq2[2].lower():
                print("Seq1 and Seq2 are identical")
            return seq1,seq2
            
def Refine_Coords(seq1,seq2,extract):
    if seq1 != "":
        print(seq1,seq2,extract)
        #print(int(seq1[0].split("-")[0])+int(seq1[3])+1,int(seq2[0].split("-")[0])-2+int(seq1[1]))
        start_after_TSD = int(seq1[0].split("-")[0])+int(seq1[3])
        end_before_TSD = int(seq2[0].split("-")[0])-1+int(seq1[1])
        RM_coords = [extract[0],start_after_TSD,end_before_TSD]
        TSDs = [extract[0],int(seq1[0].split("-")[0])+int(seq1[1])-1,int(seq1[0].split("-")[0])+int(seq1[3])-1,int(seq2[0].split("-")[0])+int(seq2[1])-1,int(seq2[0].split("-")[0])+int(seq2[3])-1]
        #print(TSDs)
    else:
        RM_coords = extract
        TSDs = []    
    return(RM_coords,TSDs)
# def ruun_water(flank, coords, extraction_path,faidx_output,ref):
        # #process the alignment
        # if seq1[2] == faidx_output[0].rstrip() or seq2[2] == faidx_output[0].rstrip():
            # print("Smith-Waterman TSD agrees with Age. Moving to next alignment")
            # same=True
            # #match +=1
            # return(seq1,seq2,between,same)
            # print("\n\n\n")
        # else:
            # print("No overlap detected. Same is false",seq1[2],seq2[2],faidx_output[0].rstrip())
            # same = False
            
        # #remove likely inaccurate alignments
        # #mismatch+=1

           
            
            # #More information needed for testing. Align the upstream and the downstream to reference at this locus
            # command = f"samtools faidx {ref} {coords[0]}:{int(coords[1])-(2*flank-1)}-{int(coords[2])+flank} > Ref_Extract.fa"
            # subprocess.run(command,shell=True)
            # cmd = f"water Ref_Extract.fa upstream_and_5_start.fa -gapopen {gap_open} --gapextend {gap_extend} -datafile Custom_Matrix -outfile Ref_and_Up.txt"
            # #print(cmd)
            # subprocess.run(cmd, shell=True)
            # upseq1,upseq2,upbetween = process_file("Ref_and_Up.txt")
            
            
            # cmd = f"water Ref_Extract.fa downstream_and_3_end.fa -gapopen {gap_open} --gapextend {gap_extend} -datafile Custom_Matrix -outfile Ref_and_Down.txt"
            # #print(cmd)
            # subprocess.run(cmd, shell=True)
            # downseq1,downseq2,downbetween = process_file("Ref_and_Down.txt")
            # print(upseq1,upbetween,upseq2)
            # print(downseq1,downbetween,downseq2)
            
            # fail = False
            # if abs(int(upseq2[3])-int(seq1[3])) > dist:
                # print(f"Locus rejected. The end of the putative 5' TSD ({seq1[3]}) is more than {dist} bp from the end of the alignment between the start of the insertion and the reference {int(upseq2[3])}.")
                # fail=True
                # #sys.exit()
            # if abs(int(downseq2[1])-int(seq2[1])) > dist:
                # print(f"Locus rejected. The start of the putative 3' TSD ({seq2[1]}) is more than {dist} bp from the end of the alignment between the end of the insertion and the reference {int(downseq2[1])}.")
                # fail=True
                # #sys.exit()
            # if fail == True:    
                # seq1[2] = "alignment_fail"
                # seq2[2] = "alignment_fail"
                # os.remove("Ref_and_Down.txt")
                # os.remove("Ref_and_Up.txt")
        # else:
            # print(f"Seq1 and Seq2 start and end within {dist} of the insertion site")
        
        # re_run = False
        # subprocess.run("rm water.txt",shell=True)
        # return(seq1,seq2,between,same)
        
        
        
def Repeat_Positions(repeat_masker_output, RM_position):
    #load in the relevant portions of the out file.
    if repeat_masker_output.endswith(".bed"):
        command = f"cat {repeat_masker_output} | grep -w '^{RM_position[0]}'"
        print(command)
        proc = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
        rm_out,err = proc.communicate()
        rm_out = rm_out.decode()
        rm_out = rm_out.split("\n")
        #print(rm_out)
        
        start_index = 1
        end_index = 2
        offset = 1
        ME_Type = 4
        O = 8
        
    else:
        command = f"cat {repeat_masker_output}"
        proc = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
        rm_out,err = proc.communicate()
        rm_out = rm_out.decode()
        rm_out = rm_out.split("\n")
        #print(rm_out)
        
        start_index = 5
        end_index = 6
        offset = 0
        ME_Type = 10
        O = 8
    
    print(f"RM for Repeat_Positions is {repeat_masker_output}")
    #new method of determining % L1 content.
    lower_bound = RM_position[1]
    upper_bound = RM_position[2]
    
    L1s = []
    #with open(repeat_masker_output, 'rt') as file:
    for line in rm_out:
        if line.startswith("There were no repetitive sequences detected"):
            break
        if line.startswith("   SW") or line.startswith("score") or line == "\n" or line == "":
            continue
        line = line.split()
        #Define relevant parameters
        #print(line)
        RM_start_in_seq = int(line[start_index])+ offset
        RM_end_in_seq = int(line[end_index])
        Element_Class = line[ME_Type]
        Element_Orientation = line[O]
        
        #create list containing the location of the sequence relative to the LINE-1
        if repeat_masker_output.endswith(".bed") and Element_Orientation == "+":
            orientation = "Forward"
            Position_within_Repeat = [int(line[9].split(',')[0]),int(line[9].split(',')[1])]
        elif repeat_masker_output.endswith(".bed") and Element_Orientation == "C":
            Position_within_Repeat = [int(line[9].split(',')[2]),int(line[9].split(',')[1])]
            orientation = "Reverse"
        elif Element_Orientation == "+":
            orientation = "Forward"
            Position_within_Repeat = [int(line[11]),int(line[12])]
        elif Element_Orientation == "C":
            Position_within_Repeat = [int(line[13]),int(line[12])]
            orientation = "Reverse"
        line_start = 0
        line_end = 0

        if Element_Class == "LINE/L1":
            #print(f"Lower and Upper Bounds of the detected insertion are {lower_bound} and {upper_bound}")
            #print(f"Start and end of detected repeat are {RM_start_in_seq} {RM_end_in_seq}")
            #print(line)
            #if RM is outside of the insertion
            if RM_end_in_seq < int(lower_bound):
                #print("L1 outide insertion")
                continue
            elif RM_start_in_seq > int(upper_bound):
                #print("L1 outide insertion")
                break
            #if RM extends beyond both ends of the insertion    
            elif RM_start_in_seq < int(lower_bound) and RM_end_in_seq > int(upper_bound):
                #print(line)
                line_start = int(lower_bound)
                line_end = int(upper_bound)
                #print("L1 extends beyond insertion boundaries")
                
            #if the RM covers one of the boundaries of the insertion
            elif RM_start_in_seq <= int(lower_bound) and RM_end_in_seq <= int(upper_bound):
                #print(line)
                line_start = int(lower_bound)
                line_end = RM_end_in_seq
                #print("L1 overlaps one boundary")
            elif RM_start_in_seq >= int(lower_bound) and RM_start_in_seq <= int(upper_bound) and RM_end_in_seq >= int(upper_bound):
                #print(line)
                line_start = RM_start_in_seq
                line_end = int(upper_bound)
                #print("L1 overlaps one boundary")
            
            #if the RM is entirely within the insertion
            elif RM_start_in_seq > int(lower_bound) and RM_start_in_seq < int(upper_bound) and RM_end_in_seq < int(upper_bound):
                #print(line)
                line_start = RM_start_in_seq
                line_end = RM_end_in_seq
                #print("L1 is within insertion boundary")
            else:
                print("Erroneous Read Detected")
                print(line)
                sys.exit("Erroneous RM Detected")
            #print("Appending L1")
            #print(line)
            L1s.append([line_start, line_end, orientation, Position_within_Repeat])
                
    #print(initial_range)
    #print(L1s)
                
                
    #combine nearby elements.
    diff = 50
    exit_loop = False
    if len(L1s) == 0:
        print("Length of L1s list is 0. Terminating")
        sys.exit()
    while not exit_loop:
        #print('hi',len(L1s))
        merged_L1 = []
        for i in range(len(L1s)):
            #print('hi')
            if i+1 == len(L1s):
                exit_loop = True
                break            
            #sys.exit()
            if L1s[i][2] == L1s[i+1][2]:
                #print(f"Orientations are the same",L1s[i][2],L1s[i+1][2])
                #sys.exit()
                if L1s[i][2] == "Forward" and abs(L1s[i+1][3][0] - L1s[i][3][1]) < diff:
                    #print(f"L1s are within range")
                    merged_L1 = [L1s[i][0],L1s[i+1][1],L1s[i][2],[L1s[i][3][0],L1s[i+1][3][1]]]
                    #print(merged_L1)
                elif L1s[i][2] == "Reverse" and abs(L1s[i][3][1] - L1s[i+1][3][0]) < diff:
                    #print(f"L1s are within range",L1s[i][3][1],L1s[i+1][3][0],diff)
                    merged_L1 = [L1s[i][0],L1s[i+1][1],L1s[i][2],[L1s[i][3][0],L1s[i+1][3][1]]]
                    #print(merged_L1)
                #else:
                    #print(f"L1s are close in contig distance, but do not overlap based on LINE-1 coordinates. The distance between the elements is larger than {diff}")

            else:
                #print(f"Orientations are different",L1s[i][2],L1s[i+1][2])
                #sys.exit()
                if max(int(L1s[i][3][0]),int(L1s[i][3][1])) > max(int(L1s[i+1][3][0]),int(L1s[i+1][3][1])):
                    #sys.exit()
                    if L1s[i][2] == "Forward":
                            #sys.exit()
                            continue
                    if abs(L1s[i][3][1] - L1s[i+1][3][1]) > diff:
                        #print("L1s are too far apart. Moving on")
                        continue
                    else:
                        #print("L1s are within range")
                        merged_L1 = [L1s[i][0],L1s[i+1][1],L1s[i][2],[L1s[i][3][0],L1s[i+1][3][0]]]
                        #print(merged_L1)
                        #sys.exit()
                    #sys.exit()
                else:
                    #diff2 = abs(max(L1s[i][0],L1s[i][1]) -  min(L1s[i+1][0],L1s[i+1][1]))
                    if L1s[i][2] == "Forward":
                        #sys.exit()
                        continue
                    if abs(L1s[i][3][0] - L1s[i+1][3][0]) > diff:
                        #print("L1s are too far apart. Moving on")
                        continue
                    else:
                        #print("L1s are within range")
                        merged_L1 = [L1s[i][0],L1s[i+1][1],L1s[i][2],[L1s[i][3][1],L1s[i+1][3][1]]]
                        #print(merged_L1)
                        #sys.exit()
            if merged_L1 != []:
                L1s[i] = merged_L1
                L1s.pop(i+1)
                break
                            #print(merged_L1)
    #print(L1s)
    #return L1s
    #sys.exit()
    
    #sys.exit()
    #determine L1 content between TSDs
    #print(f"L1s is {L1s}")
    if len(L1s) > 1:
        L1_pos = "["
        L1_content = 0
        orientations = []
        for i in range(len(L1s)):
            item = L1s[i]
            if i > 0:
                if item[0] <= L1s[i-1][1]:
                    item[0] = L1s[i-1][1] + 1
            orientations.append(item[2])
            L1_pos = L1_pos + f"[{item[3][0]},{item[3][1]}]"
            #print(f"{L1_content} + ({item[1]}-{item[0]})+1")
            L1_content = L1_content + (item[1]-item[0])+1
        L1_pos = L1_pos + "]"
        print(L1_pos)
        transduction = ["N/A","N/A"]
        if len(set(orientations)) == 1:
            orientation = orientations[0]
        else:
            orientation = "Multiple_L1_Orientations_detected"
        if orientation == "Multiple_L1_Orientations_detected":
            transduction = ["N/A","N/A"]
        elif orientation == "Forward":
            if int(L1s[-1][1]) - int(upper_bound) == 0:
                transduction = ["N/A","N/A"]
            else:
                transduction = [L1s[-1][1],upper_bound]
            #print(transduction)
        else:
            if int(lower_bound) - int(L1s[0][0]) == 0:
                transduction = ["N/A","N/A"]
            else:
                transduction = [lower_bound,L1s[0][0]]
            #print(transduction)

    else:
        if L1s == []:
            L1_content = "No_L1_detected"
            transduction = ["N/A","N/A"]
            orientation = "No_L1_detected"
            L1_pos = "No_L1_detected"
            #print('goop')
        else:
            orientation = L1s[0][2]
            L1_content = (L1s[0][1]-L1s[0][0])+1
            L1_pos = f"[{L1s[0][3][0]},{L1s[0][3][1]}]"
            print(L1_content)
            if L1_content == 1:
                L1_content = 0
                transduction = ["N/A","N/A"]
            else:
                #putative 3' transduction detection
                if orientation == "Forward":
                    if int(L1s[0][1]) - int(upper_bound) == 0:
                        transduction = ["N/A","N/A"]
                    else:
                        transduction = [L1s[0][1],upper_bound]
                    #print(transduction)
                else:
                    if int(lower_bound) - int(L1s[0][0]) == 0:
                        transduction = ["N/A","N/A"]
                    else:
                        transduction = [lower_bound,L1s[0][0]]
                    #print(transduction)
                    
    print(L1s)
    #if len(L1s) > 1:
    #sys.exit()
    return L1_content, transduction, orientation, L1_pos, L1s
def extract_for_Poly_A(transduction,extraction_path,extracted_seq,extra_seq,orientation,internal,interior_dist,exterior_cutoff=0):
    """
    Creates the samtools faidx command which will extract sequence for poly(A) detection.
        
    Args:
        transduction: Currently not used, will be needed in a future version for detecting non-reference retroelements
        extraction_path: (str) path to the fasta containing the sequence being extracted
        extracted_seq: (list) containing the sequence name, start base, and end base of the locus of interest. 1 based 
        extra_seq: (int) the number of base pairs that are extracted beyond the locus of interest 
        orientation: (str) Forward or Reverse, the orientation of the retroelement in reference to the chromosome or contig
        internal: (bool) Determines if Poly(A)s can be detected within the locus of interest. Should be set to False for retrogene detection. Later, True will be added for retrogene hallmark detection
        interior_dist: (int) The number of base pairs within the locus of interest that will be used to search for Poly(A)s. Can be used when search_internal is False
        exterior_cutoff: (int) The minimum/maximum (depending on orientation) distance that the poly(A) can be expanded into. This prevents expansion into TSDs. Default = 0
        
    Returns:
        faidx_cmd: (str) the samtools faidx command that will be used to extract sequence. If the TSD is directly adjacent to the locus of interest, no poly(A) will be found and faidx_cmd is equal to "No Poly A"
        start_in_contig: (int) The distance into the contig/chromosome of the first base extracted. 1 based.
    """
    
    if internal == True:
        if transduction[0] != "N/A":
            if orientation == "Reverse":
                #extract from end of TSD to the LINE-1
                faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{transduction[0]}-{min(int(transduction[1])+extra_seq,extracted_seq[2])}"
                start_in_contig = transduction[0]
            elif orientation == "Forward":
                faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{max(int(transduction[0])-extra_seq,extracted_seq[1])}-{int(transduction[1])}"
                start_in_contig = max(int(transduction[0])-extra_seq,extracted_seq[1])
        else:
            #don't forget to filter to not go earlier than the start of the element.
            if orientation == "Forward":
                faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{max(extracted_seq[2]-extra_seq,extracted_seq[1])}-{extracted_seq[2]}"
                start_in_contig = max(extracted_seq[2]-extra_seq,extracted_seq[1])
            elif orientation == "Reverse":
                faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{extracted_seq[1]}-{min(extracted_seq[1]+extra_seq,extracted_seq[2])}"
                start_in_contig = extracted_seq[1]
    elif internal == False:
        if exterior_cutoff == 0 and orientation == "Forward":
            exterior_cutoff = 100000000000
        if orientation == "Forward":
            faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{extracted_seq[2]-(interior_dist-1)}-{min(extracted_seq[2]+extra_seq,exterior_cutoff)}"
            extract_length = min(extracted_seq[2]+extra_seq,exterior_cutoff)-(extracted_seq[2]-(interior_dist-1))
            start_in_contig = extracted_seq[2]-(interior_dist-1)
        elif orientation == "Reverse":
            faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{max(extracted_seq[1]-extra_seq,exterior_cutoff)}-{extracted_seq[1]+(interior_dist-1)}"
            extract_length = (extracted_seq[1]+(interior_dist-1)) - max(extracted_seq[1]-extra_seq,exterior_cutoff)
            start_in_contig = max(extracted_seq[1]-extra_seq,exterior_cutoff)
        
        #occurs if the TSD starts at the first base of sequence. Thus no Poly(A) can be detected.
        if extract_length == -1:
            faidx_cmd = "No Poly A"
        elif extract_length < -1:
            sys.exit("Invalid length for poly(A) discovery")
    else:
        sys.exit()
    #print(start_in_contig)
    return faidx_cmd, start_in_contig
    
def Detect_Poly_As(transduction,extraction_path,extracted_seq,extra_seq,max_dist,orientation,internal,chrom_length,exterior_cutoff=0):
    """
    Detects Poly(A) tails located within extra_seq base pairs of a locus of interest. 
        
    Args:
        transduction: Currently not used, will be needed in a future version for detecting non-reference retroelements
        extraction_path: (str) path to the fasta containing the sequence being extracted
        extracted_seq: (list) containing the sequence name, start base, and end base of the locus of interest. 1 based 
        extra_seq: (int) the number of base pairs that are extracted beyond the locus of interest 
        max_dist: (int) detected Poly(A)s must be within dist bp of the start of the locus to be considered a valid poly(A)
        orientation: (str) Forward or Reverse, the orientation of the retroelement in reference to the chromosome or contig
        internal: (bool) Determines if Poly(A)s can be detected within the locus of interest. Should be set to False for retrogene detection. Later, True will be added for retrogene hallmark detection
        chrom_length: the length of the chromosome or contig in which the locus of interest resides
        exterior_cutoff: (int) The minimum/maximum (depending on orientation) distance that the poly(A) can be expanded into. This prevents expansion into TSDs. Default = 0
        
    Returns:
        faidx_cmd: (str) the samtools faidx command that will be used to extract sequence. If the TSD is directly adjacent to the locus of interest, no poly(A) will be found and faidx_cmd is equal to "No Poly A"
        start_in_contig: (int) The distance into the contig/chromosome of the first base extracted. 1 based.
    """
    re_run = True
    if orientation != "Forward" and orientation != "Reverse":
        print("Multiple LINE-1s in different orientations. Cannot yet resolve poly(A) or 3' transduction")
        poly_a = []
        
    else:   
        while re_run == True:
            interior_dist = 0
            re_run = False
            

            faidx_cmd, start_in_contig = extract_for_Poly_A(transduction,extraction_path,extracted_seq,extra_seq,orientation,internal,interior_dist,exterior_cutoff)
            print(faidx_cmd)
            trailing_sequence = ""
            
            if faidx_cmd == "No Poly A":
                poly_a = []
                break
                
            proc = subprocess.Popen(faidx_cmd,shell=True,stdout=subprocess.PIPE)
            faidx,err = proc.communicate()
            faidx = faidx.decode()
            
            if orientation == "Reverse":
                A_or_T = "T"
            elif orientation == "Forward":
                A_or_T = "A"
            poly_a = []
            poly_a_in_progress = ""
            start = 0
            end = None
            min_poly = 5
            for line in faidx.split("\n")[1::]:
                trailing_sequence = trailing_sequence + (line.rstrip())
            trailing_sequence = trailing_sequence.upper()
            for i in range(len(trailing_sequence)):
                if trailing_sequence[i] == A_or_T and poly_a_in_progress == "":
                    poly_a_in_progress = poly_a_in_progress + (trailing_sequence[i])
                    start = i + start_in_contig
                elif trailing_sequence[i] == A_or_T and not i == len(trailing_sequence) -1:
                    poly_a_in_progress  = poly_a_in_progress + (trailing_sequence[i])
                elif trailing_sequence[i] != A_or_T and poly_a_in_progress != "" or i == len(trailing_sequence)-1:
                    upcoming_poly = 0
                    end_of_extract = False
                    for j in range(i+1,i+6):
                        if j >= len(trailing_sequence):
                            end_of_extract = True
                            break
                        if trailing_sequence[j] == A_or_T:
                            upcoming_poly +=1
                    if upcoming_poly >= 4 or end_of_extract == True and i+1 != len(trailing_sequence):
                        poly_a_in_progress = poly_a_in_progress + trailing_sequence[i]
                    else:
                        #trim the ends of the poly if they have mutations in them.
                        if i+1 == len(trailing_sequence) and trailing_sequence[i] == A_or_T:
                            poly_a_in_progress = poly_a_in_progress + trailing_sequence[i]
                        
                        if not len(poly_a_in_progress) < min_poly:
                            end_poly = 3
                            trim_length_start = 0
                            while poly_a_in_progress[0:end_poly] != end_poly*A_or_T:
                                if len(poly_a_in_progress) <= min_poly:
                                    poly_a_in_progress = ""
                                    break
                                else:
                                    poly_a_in_progress = poly_a_in_progress[1::]
                                    start+=1
                            while poly_a_in_progress[-1:-1*end_poly-1:-1] != end_poly*A_or_T:
                                if len(poly_a_in_progress) <= min_poly:
                                    poly_a_in_progress = ""
                                    break
                                else:
                                    poly_a_in_progress = poly_a_in_progress[0:-1]
                            if len(poly_a_in_progress) > min_poly:
                                if orientation == "Forward":
                                    #start = start - extra_seq
                                    end = start + len(poly_a_in_progress)
                                    if (extracted_seq[2] + extra_seq) - (end-1) < 5 and exterior_cutoff != 0 and exterior_cutoff != "0":
                                        if extracted_seq[2] + extra_seq >= chrom_length:
                                            re_run = False
                                            print("Re-run not happening. At end of chromosome")
                                        re_run = True
                                        extra_seq +=5
                                        print(poly_a_in_progress)
                                        print("Re-run happening, Forward")
                                elif orientation == "Reverse":
                                    end = start + len(poly_a_in_progress)
                                    if abs((extracted_seq[1] - extra_seq) - (start)) < 5 and exterior_cutoff != 0 and exterior_cutoff != "0":
                                        if extracted_seq[1] - extra_seq > 1:
                                            re_run = True
                                            extra_seq +=5
                                            print(poly_a_in_progress)
                                            print("Re-run happening, Reverse")
                                        else:
                                            print("Re-run not happening. Too close to start of chrom")
                                poly_a.append([start,poly_a_in_progress,end-1])                               
                        poly_a_in_progress = ""
    print(poly_a)
    final_polys = []
    
    if internal == False:
        terminated_early = False
        if len(poly_a) > 0:
            for poly in poly_a:
                if orientation == "Reverse":
                    if poly[2] < extracted_seq[1] - max_dist:
                        print("Out of range")
                    else:
                        final_polys.append(poly)
                        if (extracted_seq[1] - extra_seq) + poly[0] < 5:
                            terminated_early = True
                elif orientation == "Forward":
                    if poly[0] > extracted_seq[2] + max_dist:
                        print("Out of range")
                    else:
                        final_polys.append(poly)
                        if (extracted_seq[2] + extra_seq) - poly[2] < 5:
                            terminated_early = True
                else:
                    sys.exit()
            print(final_polys)
    return final_polys, terminated_early
def transduction_detection(poly_a_list):
    if len(poly_a_list) < 2:
        return ["",""]
    else:
        Transduction = [poly_a_list[0][2]+1,poly_a_list[-1][0]-1]
        return Transduction