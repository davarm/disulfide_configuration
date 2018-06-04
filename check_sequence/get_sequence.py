import pandas as pd
from Bio import SeqIO

df = pd.read_csv('./../pals.csv')
pdb_list = df['PDB'].tolist()
pdb_list = list(set(pdb_list))
df['PDB'] = df['PDB'].astype(str)

fasta_file = "fasta.txt"
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

seq_dict =  {}
for record in fasta_sequences:
	#print record.id[0:6]
	seq_dict[record.id[0:6]] = record.seq 


for index,row in df.iterrows():
	pdb = row['PDB']
	pdb = pdb.upper()
	chain = row['chain1']

	#print cys1
	try:
		cys1 = row['Cys1']
		cys2 = row['Cys2']
		sequence =  seq_dict[pdb+':'+chain]
		print cys1,len(sequence), sequence[cys1-1]
		# print

	except KeyError:
		print pdb
	except IndexError:
		continue