import os
import subprocess
import sys
iterations = 100
generic_scripts_path = os.path.dirname(__file__)
prefix = "/home/blacksmi/links/kidd-lab/matt-projects/Stuff_For_Anthony/Identify_Hallmarks_From_Retrogenes/Permutation_Data_Revamp"
if not os.path.exists(prefix):
    os.makedirs(prefix)
os.chdir(prefix)

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
if not os.path.exists("logs"):
    os.mkdir("logs")
path_to_analysis_script = f"{generic_scripts_path}/Retrogene_Hallmarks_Detection.py"

for canine in canines.keys():
    path_to_fa = canines[canine][0]
    path_to_fai = path_to_fa + ".fai"
    path_to_retro_data=f"/home/blacksmi/links/kidd-lab/matt-projects/Stuff_For_Anthony/Identify_Hallmarks_From_Retrogenes/Scale_Up/{canine}/{canine}.blat_retrocopies_no_gaps.txt"
    gap_file = canines[canine][1]
    if not os.path.exists(canine):
        os.mkdir(canine)
    with open(canines[canine][1],'rt') as infile:
        with open(f"{canine}/{canine}.gaps.bed",'wt') as outfile:
            for line in infile:
                line = line.split()
                line[1] = str(int(line[1])-100)
                line[2] = str(int(line[2])+100)
                line = "\t".join(line)+'\n'
                outfile.write(line)
    cmd = "awk '{OFS=\"\\t\"}{print $1,$2-1,$3,$4,$5}'"
    cmd = f"{cmd} {path_to_retro_data} > {canine}/{canine}_0_based.bed"
    subprocess.run(cmd,shell=True)
    #print(cmd)
    #sys.exit()
    with open(f"{canine}/permute_and_run.txt",'wt') as outfile:
        for i in range(iterations):
            os.mkdir(f"{canine}/Iteration_{i}")
            cmd = f"cd Iteration_{i} && bedtools shuffle -i ../{canine}_0_based.bed -g {path_to_fai} -excl ../{canine}.gaps.bed -noOverlapping -chrom > Permuted_Output.bed"
            cmd = cmd + " && awk '{OFS=\"\\t\"}{print $1,$2+1,$3,$4,$5}' "
            cmd = f"{cmd} Permuted_Output.bed > Permuted_Output.txt && python {path_to_analysis_script} --locus_file Permuted_Output.txt --reference {path_to_fa}\n"
            outfile.write(cmd)
#sys.exit()
    if not os.path.exists(f"{canine}/logs"):
        os.mkdir(f"{canine}/logs")
    with open(f"{canine}/Permute_Commands_Driver.sh", 'wt') as file:
            header = f"#!/bin/bash\n\
#SBATCH --mail-user=blacksmi@umich.edu\n\
#SBATCH --mail-type=FAIL,ARRAY_TASKS\n\
#SBATCH --ntasks=1\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --ntasks-per-node=1\n\
#SBATCH --mem-per-cpu=12G\n\
#SBATCH --time=8:00:00\n\
#SBATCH --account=jmkidd0\n\
#SBATCH --partition=standard\n\
#SBATCH --output=logs/%x-%A_%a.out.log\n\
#SBATCH --export=ALL\n\
#SBATCH --array=1-{iterations}%{iterations}\n\n"
            
            command = f"{generic_scripts_path}/run-by-id-log.pl permute_and_run.txt logs/permute_and_run.log $SLURM_ARRAY_TASK_ID"
            file.write(header)
            file.write(command)
    subprocess.run(f"chmod ugo+x {canine}/Permute_Commands_Driver.sh", shell=True)
    os.chdir(canine)
    subprocess.run(f"sbatch Permute_Commands_Driver.sh",shell=True)
    os.chdir("..")
