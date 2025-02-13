CMD Gene Panel manipulation: 
1. Get genomic coordinates for every gene in the panel using biomaRt using GRCh37 --> output bed file
2. Get genomic coordinates for all exons from all gene isoforms --> output bed file
3. Get genomic coordinates for all exons from the longest isoform of the longest coding gene --> output bed file 
4. Get coding sequence of the shortest gene in panel based on transcript length --> output fasta file

VCF annotation: 
1. Determine how many variants PASS quality filters
2. Filter ALL variants using CMD gene panel, then filter to keep only PASS variants
3. Annonate filtered variants with allele frequencies from TopMed, then filter for only rare variants --> output tsv file with additional AF (allele frequencies) column
   
