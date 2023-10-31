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
        sys.exit("search_internal will be enabled in a future update. Please set search_internal == False")
         
        #to be added at a later date for non-reference retro-intersetions
    else:
        #sys.exit()
        if (int(seq1[1]) <= 5) or (flank_dist+interior_dist - int(seq2[3]) < 5):
            print("TSD may be longer than extracted flank. Extending flank and trying again")
            subprocess.run("rm water.txt",shell=True)
            return True, False
    return False, True
    
        
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
            sys.exit("search_internal will be enabled in a future update. Please set search_internal == False")
            #to be added at a later date for non-reference retro-intersetions
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
        sys.exit("internal == True will be enabled in a future update. Please set internal == False")        
        #to be added at a later date for non-reference retro-intersetions
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
        final_polys: (list) a list of lists where each sublist contains in order the start coordinate of each detected poly A, its sequence, the end coordinate of each poly A. All coordinates are 1 based. If no poly A is detected, an empty list is returned.
        terminated_early: (bool) Returns True if the final poly(A) would be extended if not for a limiting factor such as extending into the beginning/end of a chromosome or into a previously detected TSD. Primarily used for debugging. 
    """
    re_run = True
    if orientation != "Forward" and orientation != "Reverse":
        print("Multiple different orientations. Cannot resolve poly(A) tail")
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
                                if len(poly_a_in_progress) < min_poly:
                                    poly_a_in_progress = ""
                                    break
                                else:
                                    poly_a_in_progress = poly_a_in_progress[1::]
                                    start+=1
                            while poly_a_in_progress[-1:-1*end_poly-1:-1] != end_poly*A_or_T:
                                if len(poly_a_in_progress) < min_poly:
                                    poly_a_in_progress = ""
                                    break
                                else:
                                    poly_a_in_progress = poly_a_in_progress[0:-1]
                            if len(poly_a_in_progress) >= min_poly:
                                if orientation == "Forward":
                                    #start = start - extra_seq
                                    end = start + len(poly_a_in_progress)-1
                                    if (extracted_seq[2] + extra_seq) - (end) < 5 and (exterior_cutoff == 0 or exterior_cutoff == "0"):
                                        if extracted_seq[2] + extra_seq >= chrom_length:
                                            re_run = False
                                            print("Re-run not happening. At end of chromosome")
                                        re_run = True
                                        extra_seq +=5
                                        print(poly_a_in_progress)
                                        print("Re-run happening, Forward")
                                elif orientation == "Reverse":
                                    end = start + len(poly_a_in_progress)-1
                                    if abs((extracted_seq[1] - extra_seq) - (start)) <= 5 and (exterior_cutoff == 0 or exterior_cutoff == "0"):
                                        if extracted_seq[1] - extra_seq > 1:
                                            re_run = True
                                            extra_seq +=5
                                            print(poly_a_in_progress)
                                            print("Re-run happening, Reverse")
                                        else:
                                            print("Re-run not happening. Too close to start of chrom")
                                poly_a.append([start,poly_a_in_progress,end])                               
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
                        if abs(extracted_seq[1] - extra_seq) - poly[0] <= 5:
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
