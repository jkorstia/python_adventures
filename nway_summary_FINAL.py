#To get started:
#. conda/etc/profile.d/conda.sh 
#conda activate
#ipython
import pandas as pd
import datetime
import argparse
import os
import sys

#set up arguments
ap=argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help= "Output table from the nway tool")
ap.add_argument("-o", "--output", required=False, nargs="?", default="outfile", help="Desired name of output file, default='out'")
args=vars(ap.parse_args())


#read in the file "test_2"
# for the test file, this should have one row of headers, 99 rows of data, and 23 columns
df = pd.read_csv(args["input"], sep='\t', header=(0))

# write all print statements to a blank output file
sys.stdout=open(args["output"], "w")

now = datetime.datetime.now()
print("Summary Report")
print("Input file:", os.path.abspath("test_2"))
print("Report generated", now.strftime("%Y-%m-%d %H:%M"))
print("------------------------------------------------------------------------------")
nrows=len(df)
ncols=len(df.columns)
nspecies=int((ncols-3)/2)
species=df.columns[3:3+int(nspecies)].tolist()
print("Reference Genome:", species[0])
print()
print(int(nspecies), "species examined:")
for i in species:
	print(i)
print()
print("------------------------------------------------------------------------------")
print("Genotype Summary by Locus")
print()
print(nrows, "loci examined")
df2=df[species].apply(pd.Series.value_counts, axis=1)[['+', '-', 'N']].fillna(0)
#count the number of rows (insertions) that are found in all species (# '+' genotypes == # species)
mono=df2[df2['+']==nspecies].count()['+']
print()
print(mono, "loci are monomorphic and found in all species")
#count the number of rows (insertions) that are only found in the reference genome.
ref=df2[df2['+']==1].count()['+']
print(ref, "loci are unique to the reference species")
inf=nrows-mono-ref
print(inf, "loci are potentially phylogenetically informative")
print()
print("Of the remaining informative loci:") 
for i in range(0,nspecies-1):
	count=df2[(df2['+']!=1) & (df2['+']!=nspecies) & (df2['N']==i)].count()['+']
	print(count, "loci are missing data from", i, "species")
print()
print("------------------------------------------------------------------------------")
print("Summary by Species")
####################################################################
## this will generate the counts of insertions by species and create a summary dataframe table
##creates empty lists
p=[]
n=[]
m=[]
s=[]
#loops through each species to add their values for each item to their respective lists
for i in species:
	plus_counts=df[i].str.count("\+").sum()
	minus_counts=df[i].str.count("\-").sum()
	n_counts=df[i].str.count("N").sum()
	percent_successful=int(((plus_counts+minus_counts)/len(df))*100)
	p.append(str(plus_counts))
	m.append(str(minus_counts))
	n.append(str(n_counts))
	s.append(str(percent_successful))
# creates a dictionary from all of the lists you've created
d= {'+':p, '-':m, 'N':n, '%':s} 
# creates a dataframe from the dictionary and indexes the rows by species
df3=pd.DataFrame(d, index=species) 
print("+ = reference TE insertion found")
print("- = reference TE insertion NOT found")
print("N = presence/absence of reference TE insertion could not be determined")
print("% = percentage of genotypes successfully determined")
print()
print(df3)
print("------------------------------------------------------------------------------")
print("Program finished at", now.strftime("%Y-%m-%d %H:%M"))


	




	
	
	
	





	
