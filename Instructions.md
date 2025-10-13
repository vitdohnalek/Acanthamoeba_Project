# 🧬 RNA-seq Assembly and Quantification Pipeline

This document describes a typical **Trinity + Salmon** workflow for RNA-seq data assembly and quantification.  
It includes preprocessing, transcriptome assembly, ORF prediction, and transcript quantification steps.

---

## ⚙️ 1. Environment Setup

Load the required modules:

```bash
module load trinity/2.14.0.FULL-gcc-10.2.1-qz6xq53
module load samtools/1.14-gcc-10.2.1-oyuzddu
module load fastp
module load transdecoder
```

---

## 🔍 2. Quality Control of Raw Reads

```bash
fastqc FILE.fq -o fastqc_raw
```

---

## 🧠 3. Transcriptome Assembly with Trinity

Run Trinity using paired-end reads:

```bash
Trinity --seqType fq --left A.castL10_1.fq,A.castL11_1.fq --right A.castL10_2.fq,A.castL11_2.fq --CPU 10   --max_memory 50G --trimmomatic --output trinity_assembly
```

---

## 🧬 4. ORF Prediction with TransDecoder

Identify long open reading frames (ORFs) and predict coding regions:

```bash
TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta
```

---

## ✂️ 5. Read Trimming and Filtering

### Step 1: Trim with *fastp*

```bash
fastp -i BJ1.fastq -o BJ1.trimmed.fastq       -q 20 -u 30 --detect_adapter_for_pe       --length_required 36 -w 8       -h BJ1_fastp.html -j BJ1_fastp.json
```

### Step 2: Further Trim with *cutadapt*

```bash
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 36 -q 20 -o BJ1.cutadapt.p1.fastq BJ1.trimmed.fastq

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -O 3 -m 36 -o BJ1.cutadapt.p2.fastq BJ1.cutadapt.p1.fastq
```

### Step 3: Post-trimming Quality Check

```bash
fastqc FILE.fq -o fastqc_raw
```

---

## 🧩 6. Transcript Quantification with Salmon

### Load Salmon

```bash
module load salmon
```

### Generate Transcriptome FASTA

```bash
gffread -w transcriptome.fa -g your_genome.fasta your_annotation.gff3
```

---

### Option A: Index Without Decoys

```bash
salmon index -t transcriptome.fa -i salmon_index -k 31
```

---

### Option B: Index With Decoys

Create decoy list and build the index:

```bash
grep "^>" your_genome.fasta | cut -d " " -f 1 | sed 's/>//g' > decoys.txt
cat transcriptome.fa your_genome.fasta > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -i salmon_index_decoys
```

---

## 📊 7. Quantification (Single-End Reads)

Run quantification for each sample:

```bash
salmon quant -i salmon_index -l A -r BJ1.cutadapt.p2.fastq -o quant_BJ1 --validateMappings -p 10
```
