# Myrtaceae_NLR_workflow

**These are the steps involved in identifying the NBS-LRR (NLR) complement across the Myrtaceae family.**

### Use NLR-annotator https://doi.org/10.1104/pp.19.01273 as an initial approach.
### Use the hidden Markov model workflow as described below.

NB: The process is fully documented in hmm_process.script


1. hmmer/3.3 was used with the nhmmer function and DNA profile hidden markov models (EG_hmm) derived from *Eucalyptus grandis* alignments of the conserved nucleotide binding domain shared by Apaf-1, Resistance proteins and CED4 (NBARC) from coiled-coil (CC) NLR (here named nonTIR) and Toll interleukin-1 (TIR) NLR sequences. These EG_hmms were screened against the Myrtaceae_species_genome.fasta (where Myrtaceae_species refers to a species, for example Eucalyptus pauciflora). 

```
nhmmer EG_nonTIRhmm Myrtaceae_species_genome.fasta >Myrtaceae_species_nonTIRout
nhmmer EG_TIRhmm Myrtaceae_species_genome.fasta >Myrtaceae_species_TIRout
```

2. bed formatted files were created with the identified sequences (in the default hmmer inclusion list) in Myrtaceae_species_genome.fasta.

3. bedtools/2.29.2 was used to extract sequences from genome.

```
bedtools getfasta -s -fi Myrtaceae_species_genome.fasta -bed Myrtaceae_speciesTIR.bed -fo Myrtaceae_species_genome.TIR.fasta
bedtools getfasta -s -fi Myrtaceae_species_genome.fasta -bed Myrtaceae_speciesnonTIR.bed -fo Myrtaceae_species_genome.nonTIR.fasta
```

4. We took these outputs and created a single combined Myrtaceae_species_NBARC.fasta (from the Myrtaceae_species_genome.TIR.fasta and Myrtaceae_species_genome.nonTIR.fasta), removed duplicated sequences (using seqkit rmdup -s) and made species-specific hidden markov models in accordance with the HMMER user manual. That is, multiple sequence alignments made and output in stockholm format. These were used as input files to build the profile species-specific HMM.

```
hmmbuild -nucleic Myrtaceae_species_NBARC.hmm clustalo-Myrtaceae_species_NBARC.stockholm
```

5. Step 1,2 and 3 from the species-specific HMM extracted stranded NBARC sequences from the genome. NB. Headers for the extracted fasta have (+) for forward strand and (-) for reverse for the location on the scaffold specified.  

6. Nucleotide NBARC sequences were 6-frame translated using an adapted Bioseq script and the longest ORF frame output to fasta for downstream analysis. 

7. Headers from the NBARC_aa.fasta file were extracted and bedfiles created. nlr-annotator output bed files were reformatted to conform with hmm bed files. Both bed files were joined and sorted. Bedtools merge was used to join overlapping or abutting sequences and the subsequent bed file was used to extract nt sequences 20kb upstream and downstream of the NBARCs from Myrtaceae_species_genome.fasta. The NBARC fasta files were adjusted to remove the brackets in headers, for example Chr01:1-40000(+) changed to Chr01:1-40000+ so that they are suitable for braker2 annotations.

8. These fasta files were run through braker2/2.1.5 to annotate predicted coding regions using NLR Refseqs (aa) from Supplemetary Data https://doi.org/10.1101/2020.07.08.193961 as hints in --epmode.

```
braker.pl --cores=15 --genome=/scratch/myrtaceae/1.data/E_grandis_NBARC_20kbflanking.fasta --prot_seq=/scratch/myrtaceae/1.data/RefPlantNLR_aa.fa --ALIGNMENT_TOOL_PATH=/usr/local/genemark-es/4.59/ProtHint/bin/ --prg=ph --epmode --species=E_grandis --gff3
```

9. Amino acid outputs (stop codon '*' removed first) were submitted to Interproscan for annotation.

```
interproscan.sh -t p -appl Pfam,COILS,Gene3D -i $wdir/1.data/$input_fasta -f tsv,gff3 -d $outdir/
```

10. Annotated NB-LRR genes identified and chromosomal locations clarified against augustus_out.gff3 using scripts for these identifiers:
PF00931 = NB-ARC domain
PF08263 = LRR
PF12799 = LRR
PF13306 = LRR
PF13855 = LRR
PF13516 = LRR
Also need to search for this Gene3D output in interproscan because Pfam does not recognise all LRR signals.
G3DSA:1.10.8.430 = Gene3D LRR

Look also for TIR domains, Coils etc.
PF01582 = TIR
The scripts to do this are: run_NBARC.sh, run_NBLRR.sh


