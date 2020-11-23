#!/bin/bash

input_f="metadata/collected_input.tab"
pangolin_f="output/collected_pangolin.csv"
nextclade_f="output/collected_nextclade.csv"
fasta_f="output/collected_genomes.fasta"


echo "collecting input ..."
echo -e "batch\tkey\tvalue\tsample\tmads" > $input_f
for i in input/*.tab;
    do cat $i | awk -v hat=$(basename $i | sed 's/.tab//') '{ print hat "\t" $0 }' >> $input_f
done



echo "collecting pangolin ..."
head -n 1 $(ls output/*/*.pangolin.csv | head -n 1) > $pangolin_f
for i in output/*/*.pangolin.csv; do
	cat $i | grep -vE "^taxon" >> $pangolin_f;
done

echo "collecting nextclade ..."
head -n 1 $(ls output/*/*.nextclade.csv | head -n 1) > $nextclade_f
for i in output/*/*.nextclade.csv; do
	cat $i | grep -vE "seqName" >> $nextclade_f;
done


echo "collecting fasta ..."
cat  output/*/*.all.fasta > $fasta_f
