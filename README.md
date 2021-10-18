# nt2aa

Nucleotide to Amino Acid Translator


Input file must be cluster count csv file with the following fields: gene,pos,base,A,C,G,T,N,
These input files contain the total read counts for each nucleotide type
Usage: 
</i> python nt2aa.py -i inputfile.csv -o outputfile.csv -s start_pos_value -t threshold_aa_proportion_display <i> 

 if the -v TRUE option is used the codon preference details are displayed in the output file  
