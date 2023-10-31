import os
import subprocess
import sys

Data_Path = "/nfs/turbo/jmkiddscr/anthony-projects/retrocopy_analysis/blat_results/"
Hallmarks_Path = "/home/blacksmi/links/kidd-lab/matt-projects/Generic_Python_Scripts_Matt/Retrogene_Hallmarks_Detection.py"
generic_scripts_path = os.path.dirname(__file__)
Working_dir = "/home/blacksmi/links/kidd-lab/matt-projects/Stuff_For_Anthony/Identify_Hallmarks_From_Retrogenes/Scale_Up"

#created subdirectories 
if not os.path.exists(Working_dir):
    os.mkdir(Working_dir)

os.chdir(Working_dir)
if not os.path.exists("logs"):
    os.mkdir("logs")
    
#create commands
canines = {"china":["/nfs/turbo/jmkidddata/genomes/China_UNSW_CanFamBas_1.2/ref/China_UNSW_CanFamBas_1.2.fa",\
"/home/blacksmi/links/kidd-lab/genomes/China_UNSW_CanFamBas_1.2/ref/China_UNSW_CanFamBas_1.2.gaps.bed"], \
"mischka":["/nfs/turbo/jmkidddata/genomes/UU_Cfam_GSD_1.0/ref/UU_Cfam_GSD_1.0.fa",\
"/home/blacksmi/links/kidd-lab/genomes/UU_Cfam_GSD_1.0/ref/UU_Cfam_GSD_1.0.gaps.bed"], \
"mCanLor1.2":["/nfs/turbo/jmkidddata/genomes/mCanLor1.2/ref/mCanLor1.2.fa",\
"/home/blacksmi/links/kidd-lab/genomes/mCanLor1.2/ref/mCanLor1.2.gaps.bed"], \
"nala":["/nfs/turbo/jmkidddata/genomes/Nala_ASM864105v3/ref/Nala_ASM864105v3.fa"\
,"/home/blacksmi/links/kidd-lab/genomes/Nala_ASM864105v3/ref/Nala_ASM864105v3.gaps.bed"], \
"sandy":["/nfs/turbo/jmkidddata/genomes/Sandy_ASM325472v2/ref/Sandy_ASM325472v2.fa"\
,"/home/blacksmi/links/kidd-lab/genomes/Sandy_ASM325472v2/ref/Sandy_ASM325472v2.gaps.bed"], \
"tasha":["/nfs/turbo/jmkidddata/genomes/Dog10K_Boxer_Tasha_1.0.KP081776.1/ref/Dog10K_Boxer_Tasha_1.0.KP081776.1.fa",\
"/home/blacksmi/links/kidd-lab/genomes/Dog10K_Boxer_Tasha_1.0.KP081776.1/ref/Dog10K_Boxer_Tasha_1.0.KP081776.1.gaps.bed"], \
"wags":["/nfs/turbo/jmkidddata/genomes/Wags_PRJNA512874/ref/Wags_PRJNA512874.fa",\
"/home/blacksmi/links/kidd-lab/genomes/Wags_PRJNA512874/ref/Wags_PRJNA512874.gaps.bed"], \
"yella":["/nfs/turbo/jmkidddata/genomes/Yella_PRJNA610232/ref/Yella_PRJNA610232.fa",\
"/home/blacksmi/links/kidd-lab/genomes/Yella_PRJNA610232/ref/Yella_PRJNA610232.gaps.bed"], \
"zoey":["/nfs/turbo/jmkidddata/genomes/zoey/assemblies/2.3/ref/zoey.2.3.fa",\
"/home/blacksmi/links/kidd-lab/genomes/zoey/assemblies/2.3/gaps/zoey.2.3.gaps.bed"]}

dist_from_gaps = 100
total = len(canines.keys())
with open("Retrogene_CMDs.txt",'wt') as outfile:
    for canine in canines.keys():
        if not os.path.exists(canine):
            os.mkdir(canine)
        with open(f"{Data_Path}/{canine}.blat_retrocopies.sorted.txt","rt") as infile:
            with open(f"{canine}/{canine}.blat_retrocopies.sorted.txt",'wt') as out_loci:
                for line in infile:
                    line = line.split()
                    line[1] = str(int(line[1])-1)
                    line = "\t".join(line)+"\n"
                    out_loci.write(line)
                
        #seperate the hits near gaps from those more distal
        overlap_cmd = f"bedtools window -w {dist_from_gaps} -u -a {canine}/{canine}.blat_retrocopies.sorted.txt -b {canines[canine][1]} > {canine}/{canine}.blat_retrocopies_near_gaps.bed"
        no_overlap_cmd = f"bedtools window -w {dist_from_gaps} -v -a {canine}/{canine}.blat_retrocopies.sorted.txt -b {canines[canine][1]} > {canine}/{canine}.blat_retrocopies_no_gaps.bed"
        subprocess.run(overlap_cmd,shell=True)
        subprocess.run(no_overlap_cmd,shell=True)
        
        with open(f"{canine}/{canine}.blat_retrocopies_no_gaps.bed","rt") as infile:
            with open(f"{canine}/{canine}.blat_retrocopies_no_gaps.txt",'wt') as out_loci:
                for line in infile:
                    line = line.split()
                    line[1] = str(int(line[1])+1)
                    line = "\t".join(line)+"\n"
                    out_loci.write(line)

        cmd = f"cd {canine} && python {Hallmarks_Path} --locus_file {canine}.blat_retrocopies_no_gaps.txt --reference {canines[canine][0]}\n"
        outfile.write(cmd)

#create driver script
with open(f"Scale_Up_Retrogene_Driver.sh", 'wt') as file:
        header = f"#!/bin/bash\n\
#SBATCH --mail-user=blacksmi@umich.edu\n\
#SBATCH --mail-type=FAIL,ARRAY_TASKS\n\
#SBATCH --ntasks=1\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --ntasks-per-node=1\n\
#SBATCH --mem-per-cpu=8G\n\
#SBATCH --time=6:00:00\n\
#SBATCH --account=jmkidd0\n\
#SBATCH --partition=standard\n\
#SBATCH --output=logs/%x-%A_%a.out.log\n\
#SBATCH --export=ALL\n\
#SBATCH --array=1-{total}%{total}\n\n"
        
        command = f"{generic_scripts_path}/run-by-id-log.pl Retrogene_CMDs.txt logs/Scale_Up_Retrogene_Driver.log $SLURM_ARRAY_TASK_ID"
        file.write(header)
        file.write(command)
#fix permissions and run.
subprocess.run(f"chmod ug+x Scale_Up_Retrogene_Driver.sh", shell=True)
subprocess.run("sbatch Scale_Up_Retrogene_Driver.sh",shell=True)
