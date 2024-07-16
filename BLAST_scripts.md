# Reciprocal Best Hit EcCYP76L11 and OcCYP76L11

EcCYP76L11 protein sequence (`EcCYP76L11.fa`) vs O. coarctata proteome (`Oco.rm.prot.fa`)

```bash
makeblastdb -in Oco.rm.prot.fa -out coarctata_proteins.db -dbtype prot -parse_seqids
blastp -db coarctata_proteins.db -query EcCYP76L11.fa -out EcCYP76L11_vs_O_coarctata.txt -outfmt 6 -evalue 1e-50
```

OcCYP76L11 vs Echinochloa crus-galli

```bash
makeblastdb -in E.crus-galli_V3_GWHBDNR00000000.Protein.faa -out E.crus-galli_proteins.db -dbtype prot -parse_seqids
blastp -db E.crus-galli_proteins.db -query OcCYP76L11.fa -out OcCYP76L11_vs_E_crus-galli.txt -outfmt 6 -evalue 1e-50
```


# Quality control nBLAST

Quality control to ensure no potential momilactone genes were missed. nBLAST using the genomic locus of each momilactone gene from O. sativa against each of the Oryza species studied



To know on the results which chromosome belongs to which species, I want to append at the beginning of each chromosome the name of the species so it's easier to identify then in the blast results

```bash
sed 's/>/&O_sativa_/g' species_database/Osativa_323_v7.0.fa > species_database/Osativa_323_v7.0_modified.fa
sed 's/>/&O_punctata_/g' species_database/Oryza_punctata.Oryza_punctata_v1.2.dna.toplevel.fa > species_database/Oryza_punctata.Oryza_punctata_v1.2.dna.toplevel_modified.fa
sed 's/>/&O_officinalis_/g' species_database/Ooffic_genome.fa > species_database/Ooffic_genome_modiefied.fa
sed 's/>/&O_eichingeri_/g' species_database/Oeich_scaffolds.fasta > species_database/Oeich_scaffolds_modified.fasta
sed 's/>/&O_rhizomatis_/g' species_database/Orhiz_scaffolds.fasta > species_database/Orhiz_scaffolds_modified.fasta
sed 's/>/&O_alta_/g' species_database/O_alta_GWHAZTO00000000.genome.fasta > species_database/O_alta_GWHAZTO00000000.genome_modified.fasta
sed 's/>/&O_australiensis_/g' species_database/Oaustr_genome.fa > species_database/Oaustr_genome_modified.fa
sed 's/>/&O_coarctata_/g' species_database/Oco-Chr_genome_rm.fa > species_database/Oco-Chr_genome_rm_modified.fa
sed 's/>/&O_brachyantha_/g' species_database/O_bra_GCF_000231095.2_ObraRS2_genomic.fna > species_database/O_bra_GCF_000231095.2_ObraRS2_genomic_modified.fna
sed 's/>/&O_granulata_/g' species_database/O_granulata_GWHAAKB00000000.genome.fasta > species_database/O_granulata_GWHAAKB00000000.genome_modified.fasta
sed 's/>/&L_perrieri_/g' species_database/L_perrieri_GCA_000325765.3.fasta > species_database/L_perrieri_GCA_000325765.3_modified.fasta
```

`species_database/` is a directory that contains all the genome assemblies (fasta file) of the Oryza species used in this work

```bash
cat species_database/* > oryza_genomes.fa

makeblastdb \
    -in oryza_genomes.fa \
    -out oryza_database \
    -dbtype nucl 
```

Run blast

`rice_momi_genes/` is a directory that contains all the files of the genomic locus in fasta format for each of the MBGC genes

```bash
for gene in rice_momi_genes/*

do

gene_filename=$(basename $gene)
blastn \
    -db blast_database \
    -query $gene \
    -out blast_results/${gene_filename}_blast_result.txt \
    -outfmt 6 \
    -evalue 1e-40

done
```
