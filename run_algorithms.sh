#!/bin/sh

projectDir="$HOME/emdproject"
# echo $projectDir
dataDir="$projectDir/datasets/ECRDB61B/margin200"
# echo $dataDir
resultsDir="$projectDir/results"
# echo $resultsDir
bgmodel="$projectDir/rsat_meme_bg_ecoli_complete_genome.txt"
# echo $bgmodel

cd $dataDir
pwd

files="*.txt"
# echo $files
for file in $files
do
    filename=${file:0:-4}
    echo "Processing $file ..."
    for i in {0..10}
    do
        meme $dataDir/$file -dna -maxsize 100000000 -bfile $bgmodel -time 3600 -p 4 -mod tcm -minsites 1 -nostatus -nmotifs 5 -text -maxw 10+$i > $resultsDir/meme/$filename\_meme$i.txt
        streme --p $dataDir/$file -bfile $bgmodel --verbosity 1  --dna --totallength 100000000 --time 14400 --w 15 --seed 10+$i --thresh 0.05 --text > $resultsDir/streme/$filename\_streme$i.txt
    done
done
echo "Done"


# MotifSampler
# We made the following adjustments to the
# default parameter values. We searched five different
# motifs of a width of fifteen. The number of repeating runs
# was set to five. The background frequency model was gen-
# erated using the intergenic region sequences of all E. coli
# genome and the third order Markov model was used. We
# used the consensus score as the statistical measure for the
# quality of the predicted motifs.

#cd /home/spyros/metaptyxiako/2nd_semester/Algorithms_in_Molecular_Biology/ALGORITHMS/MotifSampler
cd "$projectDir/datasets/ECRDB61B/"
for file in margin200/*.txt; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' "$file" | sed 's/[[:space:]]*$//' > "${file%.txt}.fasta"; done

./MotifSampler -f margin200/*.fasta -b background.bg -n 5 -r 5 -w 15 -o output_200.gff -m matrix_200.txt




# BioProspector


awk '{if(NR==1) print ">ecoli\n"$0; else print}' ecoli.genome > ecoli.fasta
or 
awk '/^>/ {if (seq != "") {print seq}; printf("%s ",$0); seq=""; next;} { seq = seq $0 } END { print seq }' ecoli.txt > ecoli.fasta
correct
sed -n '1p; 2,$ {s/\s//g;H}; ${x;s/\n//g;p}' ecoli.genome > ecoli.fasta

# check the lines
wc -l myfile.txt

# dhmioyrgei txt me 2 grammes
sed -n '1p; 2,$ {H; $ {g;s/\n//g;p;q}}' yourfile.txt > newfile.txt

for file in /home/spyros/metaptyxiako/2nd_semester/Algorithms_in_Molecular_Biology/ALGORITHMS/Bioprospector/files/file_*.txt; do  
         filename=$(basename "$file");   
         sed -n '1p; 2,$ {H; $ {g;s/\n//g;p;q}}' "$filename" > "${filename%.txt}.fasta";  
         echo "${filename%.txt}.fasta";
done

# To run bioprospector

# Loop through all the txt files in margin200 directory
for file in margin200/*.txt; do

  # Get the filename without the extension
  filename=$(basename -- "$file")
  filename="${filename%.*}"

  # Loop through 10 iterations
  for i in {1..10}; do
    output="${filename}_${i}.txt"
    ./BioProspector.linux -i "$file" -f background.txt -W 15 -o "$output"
  done

done

# Combine all the results into a single file
mv *.txt results_bioprospector/


# To run MDscan

# Loop through all the txt files in margin200 directory
for file in margin200/*.txt; do

  # Get the filename without the extension
  filename=$(basename -- "$file")
  filename="${filename%.*}"

  # Loop through 10 iterations
  for i in {1..10}; do
    output="${filename}_${i}.txt"
    ./MDscan.linux -i "$file" -f background.txt -w 15 -s -o "$output"
  done

done

# Combine all the results into a single file
mv *.txt results_md/
