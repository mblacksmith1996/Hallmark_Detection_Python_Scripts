import sys
import os
import subprocess

def process_file(input_file):
    """
    Extract the query, reference, and between sequences from a smith-waterman alignment
    
    Paramters:
        input_file (str): Path to file from which the alignment is being extracted
    
    """
    with open(input_file) as infile:
                file_contents = infile.readlines()[32:-3:]
                #sys.exit()
                i = 0
                seq1 = ["","","",""]
                #between = ""
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
                    #elif i%4 == 1:
                    #    between = between + content[0]
                    elif i%4 == 2:
                        if seq2[0] == "":
                            seq2[0] = (content[0])
                            seq2[1] = (content[1])
                        seq2[2] = seq2[2] + content[2]
                        if i+2 == len(file_contents):
                            seq2[3] = (content[3])
                    i+=1
    return(seq1,seq2)#,between)

def investigate_TSD_validity(flank_dist,dist,seq1,seq2,length,orientation,post_poly_A):
    #determine if the alignments start or end in the correct place
    alignment_fail = []
    if post_poly_A == False:
        for seq in [seq1,seq2]:
            if abs(flank_dist - int(seq[1])) > dist and abs(flank_dist - int(seq[3])) > dist:
                if seq == seq1:
                    alignment_fail.append("seq1")
                else:
                    alignment_fail.append("seq2")
    else:
        if orientation == "Forward":
            #print("Forward orientation")
            if int(seq2[1]) > dist:
                alignment_fail.append(seq2)
            if abs(flank_dist - int(seq1[1])) > dist and abs(flank_dist - int(seq1[3])) > dist:   
                alignment_fail.append(seq1)
            #print(alignment_fail)
        elif orientation == "Reverse":
            print(seq1,seq2)
            if flank_dist - int(seq1[3]) > dist:
                alignment_fail.append(seq2)
            if abs(flank_dist - int(seq2[1])) > dist and abs(flank_dist - int(seq2[3])) > dist:   
                alignment_fail.append(seq1)
            #sys.exit()
    print(alignment_fail)
    if len(alignment_fail) >= 1:
        print(f"At least 1 alignment failed to be withing {dist}bp of the start/end of the inserted sequence")
        print(seq1,seq2)
        return False, False
        
    else:
        if post_poly_A == False:
            if (int(seq1[1]) < 5) or (int(seq2[1]) < 5) or (flank_dist*2 - int(seq1[3]) < 5) or (flank_dist*2 - int(seq2[3]) < 5):
                #print(length)
                if flank_dist > length/2:
                    False,False
                    #sys.exit()
                #print(seq1,seq2)
                print("TSD may be longer than extracted flank. Doubling flank and trying again")
                subprocess.run("rm water.txt",shell=True)
                return True, False
        elif orientation == "Forward":
            if (int(seq1[1]) < 5) or (flank_dist*2 - int(seq1[3]) < 5) or (flank_dist*2 - int(seq2[3]) < 5):
                if flank_dist > length/2:
                    False,False
                #print("TSD may be longer than extracted flank. Doubling flank and trying again")
                subprocess.run("rm water.txt",shell=True)
                return True, False
        elif orientation == "Reverse":
            #sys.exit()
            if (int(seq1[1]) < 5) or (int(seq2[1]) < 5) or (flank_dist*2 - int(seq2[3]) < 5):
                #print(length)
                if flank_dist > length/2:
                    False,False
                    #sys.exit()
                #print(seq1,seq2)
                print("TSD may be longer than extracted flank. Doubling flank and trying again")
                subprocess.run("rm water.txt",shell=True)
                return True, False
        return False, True

    #sys.exit()

def run_water(flank_dist, dist, coords, extraction_path,ref,orientation, post_poly_A):
    #Perform the Smith-Waterman alignment
    re_run = True
    while re_run:
        if post_poly_A != True:
            extract_up = f"{coords[0]}:{coords[1]-flank_dist}-{coords[1]+(flank_dist -1)}"
            extract_down = f"{coords[0]}:{coords[2]-(flank_dist-1)}-{coords[2]+flank_dist}"
        elif post_poly_A == True:
            if orientation == "Forward":
                extract_up = f"{coords[0]}:{coords[1]-flank_dist}-{coords[1]+(flank_dist -1)}"
                extract_down = f"{coords[0]}:{coords[2]+1}-{coords[2]+flank_dist}"
            elif orientation == "Reverse":
                extract_up = f"{coords[0]}:{coords[1]-flank_dist}-{coords[1]-1}"
                extract_down = f"{coords[0]}:{coords[2]-(flank_dist-1)}-{coords[2]+flank_dist}"

        command = f"samtools faidx {extraction_path} {extract_up} > upstream_and_5_start.fa ; samtools faidx {extraction_path} {extract_down} > downstream_and_3_end.fa"
        print(command)
        subprocess.run(command, shell=True)    
        gap_open = 10
        gap_extend = 10
        cmd = f"water upstream_and_5_start.fa downstream_and_3_end.fa -gapopen {gap_open} --gapextend {gap_extend} -datafile Custom_Matrix -outfile water.txt"
        #print(cmd)
        subprocess.run(cmd, shell=True)
        #sys.exit()
        seq1,seq2 = process_file("water.txt")   
        #print(seq1,seq2)
        re_run, valid_TSD = investigate_TSD_validity(flank_dist,dist,seq1,seq2,coords[2]-coords[1]+1,orientation,post_poly_A)
        if re_run == True:
            flank_dist = flank_dist*2
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
def extract_for_Poly_A(transduction,extraction_path,extracted_seq,extra_seq,orientation,internal):
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
        if orientation == "Forward":
            faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{extracted_seq[2]}-{extracted_seq[2]+extra_seq}"
            start_in_contig = extracted_seq[2]
        elif orientation == "Reverse":
            faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{extracted_seq[1]-extra_seq}-{extracted_seq[1]}"
            start_in_contig = extracted_seq[1]-extra_seq
    else:
        sys.exit()
    #print(start_in_contig)
    return faidx_cmd, start_in_contig
    
def Detect_Poly_As(transduction,extraction_path,extracted_seq,extra_seq,orientation,internal):
    if orientation != "Forward" and orientation != "Reverse":
        print("Multiple LINE-1s in different orientations. Cannot yet resolve poly(A) or 3' transduction")
        poly_a = []
        
    else:    
        faidx_cmd, start_in_contig = extract_for_Poly_A(transduction,extraction_path,extracted_seq,extra_seq,orientation,internal)
        trailing_sequence = ""
        proc = subprocess.Popen(faidx_cmd,shell=True,stdout=subprocess.PIPE)
        faidx,err = proc.communicate()
        faidx = faidx.decode()
        #print(faidx,start_in_contig)
        #sys.exit()
        
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
            #print(i,trailing_sequence[i])
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
                        #print("Reached end of extracted sequence")
                        end_of_extract = True
                        break
                    if trailing_sequence[j] == A_or_T:
                        upcoming_poly +=1
                #print(upcoming_poly)
                if upcoming_poly >= 4 or end_of_extract == True and i+1 != len(trailing_sequence):
                    poly_a_in_progress = poly_a_in_progress + trailing_sequence[i]
                else:
                #sys.exit()
                    #trim the ends of the poly if they have mutations in them.
                    #sys.exit()
                    if i+1 == len(trailing_sequence) and trailing_sequence[i] == A_or_T:
                        poly_a_in_progress = poly_a_in_progress + trailing_sequence[i]
                    
                    if not len(poly_a_in_progress) < min_poly:
                        #poly_a_in_progress = "TTTTTT" + poly_a_in_progress + "TTTTTAAT"
                        #print(poly_a_in_progress)
                        end_poly = 3
                        trim_length_start = 0
                        #trim_length_end = len(poly_a_in_progress)
                        while poly_a_in_progress[0:end_poly] != end_poly*A_or_T:
                            #print(poly_a_in_progress[0:end_poly])
                            if len(poly_a_in_progress) <= min_poly:
                                poly_a_in_progress = ""
                                #print(poly_a_in_progress)
                                break
                            else:
                                poly_a_in_progress = poly_a_in_progress[1::]
                                start+=1
                                #print("removing 1 nucleotide")
                            #print(poly_a_in_progress)
                        #print('broke')
                        while poly_a_in_progress[-1:-1*end_poly-1:-1] != end_poly*A_or_T:
                            if len(poly_a_in_progress) <= min_poly:
                                poly_a_in_progress = ""
                                break
                            else:
                                poly_a_in_progress = poly_a_in_progress[0:-1]
                                #print(poly_a_in_progress[-1:-1*end_poly-1:-1])
                            #sys.exit()
                            #print(poly_a_in_progress)
                        #print(poly_a_in_progress)
                        if len(poly_a_in_progress) > min_poly:
                            if orientation == "Forward":
                                #start = start - extra_seq
                                end = start + len(poly_a_in_progress)
                                #end = transduction[0]+trim_length_end - 100
                            elif orientation == "Reverse":
                                end = start + len(poly_a_in_progress)
                            poly_a.append([start,poly_a_in_progress,end-1])                               
                    poly_a_in_progress = ""
    #print(poly_a)
    return poly_a
def transduction_detection(poly_a_list):
    if len(poly_a_list) < 2:
        return ["",""]
    else:
        Transduction = [poly_a_list[0][2]+1,poly_a_list[-1][0]-1]
        return Transduction