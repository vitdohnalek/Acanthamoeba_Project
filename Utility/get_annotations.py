#Add protein annoations to your tsv
from Bio import SeqIO


annotations = {}
for seq_rec in SeqIO.parse("proteins_annotated.fasta", "fasta"):
	seq_ID = "_".join(seq_rec.id.split("|")[0].split("_")[:-1])
	annotation = seq_rec.description.split("|")[1]
	annotations[seq_ID] = annotation

tsv = "IDs\tAnnotation\tbaseMean\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
with open("DE_results.txt", "r") as f:
	for l in f:
		if l.startswith("TRINITY_"):
			tsv += l.split("\t")[0] + "\t"
			tsv += annotations[l.split("\t")[0]] + "\t"
			tsv += "\t".join(l.split("\t")[1:])

with open("DE_results_annotated.tsv", "w") as f:
	f.write(tsv)
