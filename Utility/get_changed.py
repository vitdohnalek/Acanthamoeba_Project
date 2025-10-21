#Write two separate fasta files with up and downregulated proteins
from Bio import SeqIO


upregulated = []
downregulated = []

with open("upregulated.list", "r") as f:
	for l in f:
		if l.startswith("TRINITY_"):
			upregulated.append(l.strip())

with open("downregulated.list", "r") as f:
	for l in f:
		if l.startswith("TRINITY_"):
			downregulated.append(l.strip())

upregulated_fasta = ""
downregulated_fasta = ""

for seq_rec in SeqIO.parse("proteins_annotated.fasta", "fasta"):
	seq_ID = "_".join(seq_rec.id.split("_")[:-1])
	function = seq_rec.description.split("|")[1]
	if  seq_ID in upregulated:
		upregulated_fasta += f">{seq_ID}|{function}\n{seq_rec.seq.strip('*')}\n"
		upregulated.remove(seq_ID)
	elif seq_ID in downregulated:
		downregulated_fasta += f">{seq_ID}|{function}\n{seq_rec.seq.strip('*')}\n"
		downregulated.remove(seq_ID)

with open("downregulated_proteins.fasta", "w") as f:
	f.write(downregulated_fasta)
with open("upregulated_proteins.fasta", "w") as f:
	f.write(upregulated_fasta)
