#!/bin/sh

#SBATCH -p keri 
#SBATCH -n 5
#SBATCH --mem=10000

# Samples:
# DC01_15_sw10L50 DC02_15_sw10L50 DC03_15_sw10L50 DC04_15_sw10L50 DC05_15_sw10L50 
# DS01_15_sw10L50 DS02_15_sw10L50 DS04_15_sw10L50
# DC04_17_sw10L50 SC01_17_sw10L50
# SC01_15_sw10L50 SC02_15_sw10L50 SC03_15_sw10L50 SC04_15_sw10L50 SC05_15_sw10L50
# SS01_15_sw10L50 SS02_15_sw10L50 SS05_15_sw10L50    

# view the bam file as a regular text file and name it "align":

for i in DC01_15_sw10L50 DC04_17_sw10L50 DS04_15_sw10L50 SC03_15_sw10L50 SS02_15_sw10L50 DC02_15_sw10L50 DC05_15_sw10L50 SC01_15_sw10L50 SC04_15_sw10L50 SS05_15_sw10L50 DC03_15_sw10L50 DS01_15_sw10L50 SC01_17_sw10L50 SC05_15_sw10L50 DC04_15_sw10L50 DS02_15_sw10L50 SC02_15_sw10L50 SS01_15_sw10L50; do ./samtools view $i.bam > $i.seqgenes.txt; done 


# Extract the contig on which each read was aligned on. I think the contig_id is in the third column (please verify with your own bam file before processing, and readjust # -f 3  if necessary). Create a list of contig_id named align1: 

for i in DC01_15_sw10L50 DC04_17_sw10L50 DS04_15_sw10L50 SC03_15_sw10L50 SS02_15_sw10L50 DC02_15_sw10L50 DC05_15_sw10L50 SC01_15_sw10L50 SC04_15_sw10L50 SS05_15_sw10L50 DC03_15_sw10L50 DS01_15_sw10L50 SC01_17_sw10L50 SC05_15_sw10L50 DC04_15_sw10L50 DS02_15_sw10L50 SC02_15_sw10L50 SS01_15_sw10L50; do cut -f 3 $i.seqgenes.txt > $i.allgenes.txt; done

# Sort the contig list and count how many time each contig_id appears in the list:

for i in DC01_15_sw10L50 DC04_17_sw10L50 DS04_15_sw10L50 SC03_15_sw10L50 SS02_15_sw10L50 DC02_15_sw10L50 DC05_15_sw10L50 SC01_15_sw10L50 SC04_15_sw10L50 SS05_15_sw10L50 DC03_15_sw10L50 DS01_15_sw10L50 SC01_17_sw10L50 SC05_15_sw10L50 DC04_15_sw10L50 DS02_15_sw10L50 SC02_15_sw10L50 SS01_15_sw10L50; do sort $i.allgenes.txt | uniq -c > $i.countgenes.txt; done

# Delete "space characters" at the beggining of each line and change the "space" between the number of contig occurrences and the contig_id for a "tab" 

for i in DC01_15_sw10L50 DC04_17_sw10L50 DS04_15_sw10L50 SC03_15_sw10L50 SS02_15_sw10L50 DC02_15_sw10L50 DC05_15_sw10L50 SC01_15_sw10L50 SC04_15_sw10L50 SS05_15_sw10L50 DC03_15_sw10L50 DS01_15_sw10L50 SC01_17_sw10L50 SC05_15_sw10L50 DC04_15_sw10L50 DS02_15_sw10L50 SC02_15_sw10L50 SS01_15_sw10L50; do awk '{ sub(/^[ \t]+/, ""); print }' $i.genesorder.txt | sed 's/ /\t/' >; done 

# There we go, you should obtain a list looking like that:

# 59 AB_000002_T.1
# 23 AB_000003_T.1
# 6 AB_000006_T.1
# 14 AB_000007_T.1

# Eliminar archivos .bam de Programas

rm *.bam

# Mover todos los archivos creados a su directorio correspondiente 

for i in DC01_15_sw10L50 DC04_17_sw10L50 DS04_15_sw10L50 SC03_15_sw10L50 SS02_15_sw10L50 DC02_15_sw10L50 DC05_15_sw10L50 SC01_15_sw10L50 SC04_15_sw10L50 SS05_15_sw10L50 DC03_15_sw10L50 DS01_15_sw10L50 SC01_17_sw10L50 SC05_15_sw10L50 DC04_15_sw10L50 DS02_15_sw10L50 SC02_15_sw10L50 SS01_15_sw10L50; do mv $i.seqgenes.txt ../../TRANSCRIPTOMICS_MAP/Count/alignment_AbP_paired_sw10_L50; done
for i in DC01_15_sw10L50 DC04_17_sw10L50 DS04_15_sw10L50 SC03_15_sw10L50 SS02_15_sw10L50 DC02_15_sw10L50 DC05_15_sw10L50 SC01_15_sw10L50 SC04_15_sw10L50 SS05_15_sw10L50 DC03_15_sw10L50 DS01_15_sw10L50 SC01_17_sw10L50 SC05_15_sw10L50 DC04_15_sw10L50 DS02_15_sw10L50 SC02_15_sw10L50 SS01_15_sw10L50; do mv $i.allgenes.txt ../../TRANSCRIPTOMICS_MAP/Count/alignment_AbP_paired_sw10_L50; done
for i in DC01_15_sw10L50 DC04_17_sw10L50 DS04_15_sw10L50 SC03_15_sw10L50 SS02_15_sw10L50 DC02_15_sw10L50 DC05_15_sw10L50 SC01_15_sw10L50 SC04_15_sw10L50 SS05_15_sw10L50 DC03_15_sw10L50 DS01_15_sw10L50 SC01_17_sw10L50 SC05_15_sw10L50 DC04_15_sw10L50 DS02_15_sw10L50 SC02_15_sw10L50 SS01_15_sw10L50; do mv $i.countgenes.txt ../../TRANSCRIPTOMICS_MAP/Count/alignment_AbP_paired_sw10_L50; done
for i in DC01_15_sw10L50 DC04_17_sw10L50 DS04_15_sw10L50 SC03_15_sw10L50 SS02_15_sw10L50 DC02_15_sw10L50 DC05_15_sw10L50 SC01_15_sw10L50 SC04_15_sw10L50 SS05_15_sw10L50 DC03_15_sw10L50 DS01_15_sw10L50 SC01_17_sw10L50 SC05_15_sw10L50 DC04_15_sw10L50 DS02_15_sw10L50 SC02_15_sw10L50 SS01_15_sw10L50; do mv $i.genesorder.txt ../../TRANSCRIPTOMICS_MAP/Count/alignment_AbP_paired_sw10_L50; done

 
