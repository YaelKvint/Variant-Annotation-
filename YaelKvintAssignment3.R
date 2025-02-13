library(biomaRt)

# Load mart 37
mart_37 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://grch37.ensembl.org")
mart_37 <- useDataset(mart_37, dataset = "hsapiens_gene_ensembl")

## Using BiomaRT to get gene coordinates for every gene in the panel (panel.bed)
# Load gene panel file 
cmd_panel <- "Congenital muscular dystrophy.tsv"
genes <- read.delim(cmd_panel, stringsAsFactors = F)

# Extract gene coordinates for all genes using their ensembl ID
genes_info <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position",
                 "ensembl_gene_id"),
  filters = c("ensembl_gene_id"), # filtering using ensembl ID
  values = genes$EnsemblId.GRch37., # ensembl IDs from the panel are used
  mart = mart_37)

# Order genes by chromosome number
genes_info_sorted <- genes_info[order(genes_info$chromosome_name, genes_info$start_position), ]

# Write sorted genes and associated info into a bed file
write.table(genes_info_sorted, 
            file = 'panel.bed',
            sep = "\t", quote = F, row.names = F, col.names = F)

#------------------------------------------------------------------------------#
## Using BiomaRT to get exon coordinates from all isoforms of every gene in the panel (panel_exons.bed)
# Get exon isoforms for all genes 
genes_exons <- getBM(
  attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "ensembl_exon_id",
                 "ensembl_gene_id"),
  filters = c("ensembl_gene_id"), # filtering using ensembl ID
  values = genes$EnsemblId.GRch37., # ensembl IDs from the panel are used
  mart = mart_37)

# Order genes by chromosome number
genes_exons_sorted <- genes_exons[order(genes_exons$chromosome_name, genes_exons$exon_chrom_start), ]

# Write exons and associated info into a bed file
write.table(genes_exons_sorted, 
            file = 'panel_exons.bed',
            sep = "\t", quote = F, row.names = F, col.names = F)

#------------------------------------------------------------------------------#
## Using BiomaRT to get all exons from the longest isoform of the longest coding gene in the panel (panel_exons.bed)
# Retrieve coding sequence lengths for all genes
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "cds_length"),
  filters = "ensembl_gene_id", # filtering using ensembl ID
  values = genes$EnsemblId.GRch37., # ensembl IDs from panel are used
  mart = mart_37
)

# Find the longest coding gene
longest_gene <- gene_info[which.max(gene_info$cds_length), "ensembl_gene_id"]

# Retrieve all gene transcripts from longest coding gene
transcripts <- getBM(
  attributes = c("ensembl_transcript_id", "cds_length"),
  filters = "ensembl_gene_id", # filtering using ensembl ID
  values = longest_gene, # specifying for the variable that contains the ensembl ID of the longest coding gene
  mart = mart_37
)

# Find the transcript with the longest coding sequence 
longest_transcript <- transcripts[which.max(transcripts$cds_length), "ensembl_transcript_id"]

# Retrieve all exons for the longest transcript
exons <- getBM(
  attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "ensembl_exon_id"),
  filters = "ensembl_transcript_id", # filtering using ensembl ID
  values = longest_transcript, # specifying for the variable that contains the ensembl ID of the longest transcript
  mart = mart_37
)

# Write exons and associated info into file 
write.table(exons, 
            file = 'panel_exons_longest_transcript.bed',
            sep = "\t", quote = F, row.names = F, col.names = F)

#------------------------------------------------------------------------------#
# Extracting the coding sequence of the shortest gene in the panel
# Retrieve transcript length and CDS information
transcript_data <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length", "cds_length"),
  filters = "ensembl_gene_id", # filtering using ensembl ID
  values = genes$EnsemblId.GRch37., # ensembl IDs from the panel are used
  mart = mart_37
)

# Filter out non-coding transcripts
coding_transcripts <- transcript_data[!is.na(transcript_data$cds_length), ]

# Find the shortest transcript with a CDS
shortest_transcript <- coding_transcripts[which.min(coding_transcripts$transcript_length), ]

# Get the CDS for the shortest transcript
cds_sequence <- getSequence(
  id = shortest_transcript$ensembl_transcript_id, # specifying ensembl ID from shortest transcript variable
  type = "ensembl_transcript_id", 
  seqType = "coding", # extracting specifically the coding sequence
  mart = mart_37
)

# Writing CDS to a fasta file 
write("> ENSG00000198947", "ENSG00000198947.fasta")
write(cds_sequence$coding, "ENSG00000198947.fasta", append = T)

#------------------------------------------------------------------------------#
## Determining how many variants PASS quality filters from downloaded vcf
library(VariantAnnotation)

# Reading in the VCF file 
vcf <- readVcf("PGPC_0001_S1.flt.subset.vcf", "hg19")

# Extracting which variants in the vcf PASS quality filters 
variants_pass <- rowRanges(vcf)[rowRanges(vcf)$FILTER == "PASS", ]
length(variants_pass) # Getting the exact number of how many variants passed 
# Number of variants passed = 42535

#------------------------------------------------------------------------------#
## Filtering for PASS variants that are also in the gene panel (panel.bed)
# Read panel.bed into R and give it column names:
panel <- read.table("panel.bed", sep = "\t", 
                    col.names = c("chr", "start", "end", "gene"))

# Transform already loaded in gene panel into a GRanges object
panel_region <- GRanges(seqnames = panel$chr,
                        ranges = IRanges(start = panel$start,
                                         end = panel$end,
                                         names = panel$gene))

# Update chromosome names to match vcf format
seqlevelsStyle(panel_region) <- "UCSC"

# Bgzip and index the vcf file
bgzip("PGPC_0001_S1.flt.subset.vcf", "PGPC_0001_S1.flt.subset.vcf.gz", 
      overwrite = T)
indexTabix("PGPC_0001_S1.flt.subset.vcf.gz", format = "vcf")

# Read index into R:
tbx <- TabixFile("PGPC_0001_S1.flt.subset.vcf.gz")

# Subset all variants in the vcf based on the panel region
vcf_subset <- readVcf(tbx, "hg19", param = panel_region)

# Filter subsetted vcf for PASS variants only
vcf_subset_pass <- vcf_subset[rowRanges(vcf_subset)$FILTER == "PASS", ]
length(vcf_subset_pass) # There are 125 PASS variants in panel genes 

# Expand the vcf so each entry has as single ALT allele
vcf_subset_expanded <- expand(vcf_subset_pass)
length(vcf_subset_expanded) # After expanding, there are 126 PASS variants in panel genes

# Save the subsetted expanded PASS variants to a TSV file
write.table(rowRanges(vcf_subset_expanded), file = "PGP_variants.tsv", 
            sep = "\t", quote = F, row.names = F)

#------------------------------------------------------------------------------#
## Annotating filtered panel variants with allele frequencies from TopMed and keeping only rare variants
library("MafDb.TOPMed.freeze5.hg19")
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")

# Create shorter variable names for the needed dataset
mafdb <- MafDb.TOPMed.freeze5.hg19 # minor allele frequency database

# Convert the expanded subset into a GRanges object for further use 
gr_variants <- rowRanges(vcf_subset_expanded)

# Get the minor allele frequency for the SNPs, as found in the TOPMed database
snp_pos_maf <- gscores(mafdb, gr_variants)

# Remove any NAs from the allele frequency column, then filter for only rare alleles (<0.01)
snp_pos_maf_noNA <- snp_pos_maf[!is.na(snp_pos_maf$AF), ]
snp_pos_maf_rare <- snp_pos_maf_noNA[snp_pos_maf_noNA$AF < 0.01, ] # Filtering for variants that have an AF of <0.0

# There are 55 variants that PASS and have an MAF < 0.01. 

# Export the rare SNPs with AF as a tsv file
write.table(as.data.frame(snp_pos_maf_rare), file = "PGP_rare_variants.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

