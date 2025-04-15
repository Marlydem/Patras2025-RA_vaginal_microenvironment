conda activate qiime2-2023.5

cd /Users/marlydmejia/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Desktop\ Documents-\ Patras\ Lab/Proposals/Manuscripts/Mine/2024_Rheumatoid_Arthritis_Vaginal

#download Greengenes2 files into parent folder so I can reach into them regardless of the study I am working on

#mkdir tables
#mkdir repseqs

#mkdir VSP80_2024updated
mkdir VSP80_2024updated/Output
mkdir VSP80_2024updated/Visualization
mkdir VSP80_2024updated/exported


#MULTIPLEXED SEQS (already done, skip to merge seqs)


gzip VSP_Human_swab_samples/2022_04_21_040122SO27F_Raw_Data-mapping/analysis-040122SO27FSam1-21-115-m64277e_220413_194908.hifi_reads.fastq
gzip VSP_Human_swab_samples/2022_04_21_040122SO27F_Raw_Data-mapping/analysis-040122SO27FSam1-21-127-m64277e_220414_181323.hifi_reads.fastq
gzip VSP_Human_swab_samples/2022_04_21_040122SO27F_Raw_Data-mapping/analysis-040122SO27FSam1-21-141-m64277e_220417_231829.hifi_reads.fastq

gzip VSP_Human_swab_samples/2022_07_15_061022CR27F_Raw_Data-mapping/m64277e_220714_152820.hifi_reads.bc1001--bc1001.bam.fastq

sed 's/,/\t/g' VSP_Human_swab_samples/2022_07_15_061022CR27F_Raw_Data-mapping/061022CR27F-mapping.csv VSP_Human_swab_samples/2022_07_15_061022CR27F_Raw_Data-mapping/061022CR27F-mapping.tsv

sed 's/,/\t/g' VSP_Human_swab_samples/2022_04_21_040122SO27F_Raw_Data-mapping/040122SO27F-mapping.csv VSP_Human_swab_samples/2022_04_21_040122SO27F_Raw_Data-mapping/040122SO27F-mapping.tsv

qiime tools import \
--type MultiplexedSingleEndBarcodeInSequence \
--input-path VSP_Human_swab_samples/2022_07_15_061022CR27F_Raw_Data-mapping/m64277e_220714_152820.hifi_reads.bc1001--bc1001.bam.fastq.gz \
--output-path multiplexed-seqs_062022.qza
#Imported VSP_Human_swab_samples/2022_07_15_061022CR27F_Raw_Data-mapping/m64277e_220714_152820.hifi_reads.bc1001--bc1001.bam.fastq.gz as MultiplexedSingleEndBarcodeInSequenceDirFmt to multiplexed-seqs_062022.qza

qiime tools import \
--type MultiplexedSingleEndBarcodeInSequence \
--input-path VSP_Human_swab_samples/2022_04_21_040122SO27F_Raw_Data-mapping/analysis-040122SO27FSam1-21-141-m64277e_220417_231829.hifi_reads.fastq.gz \
--output-path VSP01_20_MVS2023/multiplexed-seqs_032022.qza
Imported VSP_Human_swab_samples/2022_04_21_040122SO27F_Raw_Data-mapping/analysis-040122SO27FSam1-21-141-m64277e_220417_231829.hifi_reads.fastq.gz as MultiplexedSingleEndBarcodeInSequenceDirFmt to VSP01_20_MVS2023/multiplexed-seqs_032022.qza

qiime cutadapt demux-single \
--i-seqs VSP01_20_MVS2023/multiplexed-seqs_062022.qza \
--m-barcodes-file VSP_Human_swab_samples/2022_07_15_061022CR27F_Raw_Data-mapping/061022CR27F-mapping.tsv \
--m-barcodes-column BarcodeSequence \
--p-error-rate 0 \
--o-per-sample-sequences demultiplexed-seqs_MVS062022.qza \
--o-untrimmed-sequences untrimmed_MVS062022.qza \
--verbose

qiime cutadapt demux-single \
--i-seqs VSP01_20_MVS2023/multiplexed-seqs_032022.qza \
--m-barcodes-file VSP_Human_swab_samples/2022_04_21_040122SO27F_Raw_Data-mapping/040122SO27F-mapping.tsv \
--m-barcodes-column BarcodeSequence \
--p-error-rate 0 \
--o-per-sample-sequences demultiplexed-seqs_MVS032022.qza \
--o-untrimmed-sequences untrimmed_MVS032022.qza \
--verbose
#produces demultiplexed reads

qiime dada2 denoise-ccs \
--i-demultiplexed-seqs VSP01_20_MVS2023/demultiplexed-seqs_MVS062022.qza \
--p-front AGRGTTTGATCMTGGCTCAG \
--p-adapter GGGTTACCTTGTTACGACTT \
--p-max-mismatch 2 \
--p-indels FALSE \
--p-trunc-len 0 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-min-len 1000 \
--p-max-len 1600 \
--p-pooling-method independent \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 3.5 \
--p-allow-one-off FALSE \
--p-n-threads 1 \
--p-n-reads-learn 1000000 \
--p-hashed-feature-ids TRUE \
--o-table MVS062022_feature-table.qza \
--o-representative-sequences MVS062022_rep-seqs.qza \
--o-denoising-stats VSP01_20_MVS2023/Output/MVS062022_denoising-stats \
--verbose

qiime dada2 denoise-ccs \
--i-demultiplexed-seqs VSP01_20_MVS2023/demultiplexed-seqs_MVS032022.qza \
--p-front AGRGTTTGATCMTGGCTCAG \
--p-adapter GGGTTACCTTGTTACGACTT \
--p-max-mismatch 2 \
--p-indels FALSE \
--p-trunc-len 0 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-min-len 1000 \
--p-max-len 1600 \
--p-pooling-method independent \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 3.5 \
--p-allow-one-off FALSE \
--p-n-threads 1 \
--p-n-reads-learn 1000000 \
--p-hashed-feature-ids TRUE \
--o-table MVS032022_feature-table.qza \
--o-representative-sequences MVS032022_rep-seqs.qza \
--o-denoising-stats VSP01_20_MVS2023/Output/MVS032022_denoising-stats \
--verbose


#DEMULITPLEXED


qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path VSP80_2024updated/manifest-VSP80run_demultiplexed_extra.txt \
--output-path VSP80_2024updated/Single-end-demux_VSP80_demultiplexed_extra.qza \
--input-format SingleEndFastqManifestPhred33V2

qiime dada2 denoise-ccs \
--i-demultiplexed-seqs VSP80_2024updated/Single-end-demux_VSP80_demultiplexed_extra.qza \
--p-front AGRGTTTGATCMTGGCTCAG \
--p-adapter GGGTTACCTTGTTACGACTT \
--p-max-mismatch 2 \
--p-indels FALSE \
--p-trunc-len 0 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-min-len 1000 \
--p-max-len 1600 \
--p-pooling-method independent \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 3.5 \
--p-allow-one-off FALSE \
--p-n-threads 1 \
--p-n-reads-learn 1000000 \
--p-hashed-feature-ids TRUE \
--o-table VSP80_2024updated/VSP80_2024updated_feature-table_extra.qza \
--o-representative-sequences VSP80_2024updated/VSP80_2024updated_rep-seqs_extra.qza \
--o-denoising-stats VSP80_2024updated/Output/VSP80_2024updated_denoising-stats_extra \
--verbose

qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path VSP80_2024updated/manifest-VSP80run_demultiplexed.txt \
--output-path VSP80_2024updated/Single-end-demux_VSP80_demultiplexed.qza \
--input-format SingleEndFastqManifestPhred33V2

qiime dada2 denoise-ccs \
--i-demultiplexed-seqs VSP80_2024updated/Single-end-demux_VSP80_demultiplexed.qza \
--p-front AGRGTTTGATCMTGGCTCAG \
--p-adapter GGGTTACCTTGTTACGACTT \
--p-max-mismatch 2 \
--p-indels FALSE \
--p-trunc-len 0 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-min-len 1000 \
--p-max-len 1600 \
--p-pooling-method independent \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 3.5 \
--p-allow-one-off FALSE \
--p-n-threads 1 \
--p-n-reads-learn 1000000 \
--p-hashed-feature-ids TRUE \
--o-table VSP80_2024updated/VSP80_2024updated_feature-table.qza \
--o-representative-sequences VSP80_2024updated/VSP80_2024updated_rep-seqs.qza \
--o-denoising-stats VSP80_2024updated/Output/VSP80_2024updated_denoising-stats \
--verbose

qiime feature-table tabulate-seqs \
--i-data VSP80_2024updated/VSP80_2024updated_rep-seqs.qza \
--o-visualization VSP80_2024updated/Visualization/VSP80_repseqs

qiime tools export \
--input-path VSP80_2024updated/VSP80_2024updated_feature-table.qza \
--output-path VSP80_2024updated/exported/merged_VSP80_table

qiime tools export \
--input-path VSP80_2024updated/VSP80_2024updated_rep-seqs.qza \
--output-path VSP80_2024updated/exported/merged_VSP80_rep-seqs

qiime metadata tabulate \
--m-input-file VSP80_2024updated/Output/VSP80_2024updated_denoising-stats.qza \
--o-visualization VSP80_2024updated/Visualization/VSP80_2024updated_denoising-stats.qzv

qiime feature-table filter-samples \
--i-table VSP80_2024updated/VSP80_2024updated_feature-table_extra.qza \
--m-metadata-file VSP80_2024updated/manifest-VSP80run_demultiplexed_extra_filter.txt \
--o-filtered-table VSP80_2024updated/VSP80_2024updated_feature-table_extraFILTERED.qza


#MERGING ALL
qiime feature-table merge \
--i-tables VSP80_2024updated/MVS032022_feature-table.qza \
--i-tables VSP80_2024updated/MVS062022_feature-table.qza \
--i-tables VSP80_2024updated/VSP80_2024updated_feature-table.qza \
--i-tables VSP80_2024updated/VSP80_2024updated_feature-table_extraFILTERED.qza \
--o-merged-table VSP80_2024updated/FULL_VSP80_table_extra.qza
#Saved FeatureTable[Frequency] to:


qiime tools export \
--input-path VSP80_2024updated/FULL_VSP80_table_extra.qza \
--output-path VSP80_2024updated/exported/merged_VSP80_table_extra
biom convert -i VSP80_2024updated/exported/merged_VSP80_table_extra/feature-table.biom -o VSP80_2024updated/exported/merged_VSP80_table_extra/OTU_table.txt --to-tsv --header-key taxonomy

qiime feature-table merge-seqs \
--i-data VSP80_2024updated/VSP80_2024updated_rep-seqs.qza \
--i-data VSP80_2024updated/MVS032022_rep-seqs.qza \
--i-data VSP80_2024updated/MVS062022_rep-seqs.qza \
--i-data VSP80_2024updated/VSP80_2024updated_rep-seqs_extra.qza \
--o-merged-data VSP80_2024updated/FULL_VSP80_repseqs_extra
#Saved FeatureData[Sequence] to: repseqs/VSP_MVS_1.qza

qiime feature-table filter-seqs \
--i-data VSP80_2024updated/FULL_VSP80_repseqs_extra.qza \
--i-table VSP80_2024updated/FULL_VSP80_table_extra.qza \
--o-filtered-data VSP80_2024updated/FULL_VSP80_repseqs_extra_filtered

qiime feature-table tabulate-seqs \
--i-data VSP80_2024updated/FULL_VSP80_repseqs_extra.qza \
--o-visualization VSP80_2024updated/Visualization/FULL_VSP80_repseqs_extra

qiime tools export \
--input-path VSP80_2024updated/VSP80_2024updated_rep-seqs_extra.qza \
--output-path VSP80_2024updated/exported/FULL_VSP80_repseqs_extra




#start taxonomy
pip install q2-greengenes2

#non-V4 because I have FULL-length sequences
qiime greengenes2 non-v4-16s \
--i-table VSP80_2024updated/FULL_VSP80_table_extra.qza \
--i-sequences VSP80_2024updated/FULL_VSP80_repseqs_extra_filtered.qza \
--i-backbone /Users/marlydmejia/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Desktop\ Documents-\ Patras\ Lab/Proposals/Manuscripts/Mine/2024_Rheumatoid_Arthritis_Vaginal/VSP83_2024/2022.10.backbone.full-length.fna.qza \
--o-mapped-table VSP80_2024updated/Casecont1_feattable_gg2.qza \
--o-representatives VSP80_2024updated/VSP80_seqs2_gg2.qza
#Saved FeatureTable[Frequency] to:
#Saved FeatureData[Sequence] to:

qiime feature-table tabulate-seqs \
--i-data VSP80_2024updated/VSP80_seqs2_gg2.qza \
--o-visualization VSP80_2024updated/Visualization/VSP80_seqs2_gg2

qiime tools export \
--input-path VSP80_2024updated/Casecont1_feattable_gg2.qza \
--output-path VSP80_2024updated/exported/Casecont1_feattable_gg2
biom convert -i VSP80_2024updated/exported/Casecont1_feattable_gg2/feature-table.biom -o VSP80_2024updated/exported/Casecont1_feattable_gg2/OTU_table.txt --to-tsv --header-key taxonomy


qiime feature-table filter-samples \
--i-table VSP80_2024updated/Casecont1_feattable_gg2.qza \
--m-metadata-file VSP80_2024updated/samples-to-keep_VSP87MVS.txt \
--o-filtered-table VSP80_2024updated/FULL_VSP80_gg2table


qiime greengenes2 taxonomy-from-table \
--i-reference-taxonomy VSP83_2024/2022.10.taxonomy.md5.nwk.qza \
--i-table VSP80_2024updated/FULL_VSP80_gg2table.qza \
--o-classification VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza

qiime metadata tabulate \
--m-input-file VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza \
--o-visualization VSP80_2024updated/Visualization/VSP80_01_gg2.tabletaxonomymd5.qza

##EXPORTING##
qiime tools export \
--input-path VSP80_2024updated/FULL_VSP80_gg2table.qza \
--output-path VSP80_2024updated/exported/FULL_VSP80_gg2table

qiime tools export \
--input-path VSP80_2024updated/VSP80_seqs2_gg2.qza \
--output-path VSP80_2024updated/exported/VSP80_seqs2_gg2

qiime tools export \
--input-path VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza \
--output-path VSP80_2024updated/exported/VSP80_01_gg2.tabletaxonomymd5




#OVERALL peek into table to see what taxa are in blanks and need to be removed
biom convert -i VSP80_2024updated/exported/FULL_VSP80_gg2table/feature-table.biom -o VSP80_2024updated/exported/FULL_VSP80_gg2table/OTU_table.txt --to-tsv --header-key taxonomy
  #biom convert -i VSP80_2024updated/exported/FULL_VSP80_gg2table/OTU_gg2taxatable.txt -o VSP80_2024updated/exported/FULL_VSP80_gg2table/OTU_TAXAtable_json.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy

#------------------------------------------------------------------------------------------------------------------------------------------------------------
#didn't use this time, but here for reference
qiime quality-control decontam-identify \
--i-table FULL_VSP80_gg2table.qza \
--m-metadata-file Metadata_RAstudy_ControlVector.csv \
--p-method 'prevalence' \
--p-prev-control-column 'Sample_or_Control' \
--p-prev-control-indicator 'Control' \
--o-decontam-scores Decontam_scores

#WOULD WANT TO CAPTURE UNIQUE SIGNATURES IN INDIVIDUALS, SO FOUND IN ONE SAMPLE IS ACCEPTABLE

#------------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%% R studio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ######### FILTER OUT CONTAMINANTS USING DECONTAM #########

## ----loadPS----------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

" skipped this chunck since following warning came up
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
above three lines results in
'getOption("repos")' replaces Bioconductor standard repositories, see 'help("repositories",
package = "BiocManager")' for details.
Replacement repositories:
    CRAN: https://cran.rstudio.com/
Bioconductor version 3.18 (BiocManager 1.30.22), R 4.3.1 (2023-06-16)
Warning message:
package(s) not installed when version(s) same as or greater than current; use `force = TRUE` to
  re-install: 'decontam'
"

library(phyloseq); packageVersion("phyloseq")
"1.46.0"

library(ggplot2); packageVersion("ggplot2")
"3.5.1" #upgraded from 3.5.0

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")
library(decontam); packageVersion("decontam")
"1.22.0"

setwd("~")
#upload the matrix into phyloseq format
#transpose feat table and change first col title to "sample_name"
#Go in and add prefix "ID_" to all feature ID names in the taxonomy file and in the OTU tables and such IF ANY TAXONOIC IDENTIFIERS START WITH A NUMBER. This way they can be matched later for taxa removal or merging. IF THERE ARE HYPHENS AS IN WITH gg2022, CONVERT THEM TO PERIODS (will make it hard to target the periods and change them back as periods and other symbols exist)
#feature table read with feature IDs and not taxonomy names
#phyloseq converts "-" to periods, so need to replace characters in original Taxonomy_IDs file prior to import
#make sure "taxonomy" column (now last row) is deleted
#add col in metadata "Sample_or_Control" and indicate Control or True Sample

tFeatureTable <- read.csv("/Users/marlydmejia/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Desktop\ Documents-\ Patras\ Lab/Proposals/Manuscripts/Mine/2024_Rheumatoid_Arthritis_Vaginal/VSP80_2024updated/exported/FULL_VSP80_gg2table/OTU_fliptaxa_table.csv", header=TRUE, sep = ",", row.names="sample_name")

metadata_set <- read.csv("/Users/marlydmejia/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Desktop\ Documents-\ Patras\ Lab/Proposals/Manuscripts/Mine/2024_Rheumatoid_Arthritis_Vaginal/Metadata_RAstudy_ControlVector.csv", header=TRUE, sep = ",", row.names="sample_name")

taxonomy <- read.csv("/Users/marlydmejia/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Desktop\ Documents-\ Patras\ Lab/Proposals/Manuscripts/Mine/2024_Rheumatoid_Arthritis_Vaginal/VSP80_2024updated/exported/VSP80_01_gg2.tabletaxonomymd5/taxonomy_flip.csv", header=TRUE, sep = ",", row.names="Feature_ID")

#create directory "Decontam_process"
dir.create("/Users/marlydmejia/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Desktop\ Documents-\ Patras\ Lab/Proposals/Manuscripts/Mine/2024_Rheumatoid_Arthritis_Vaginal/VSP80_2024updated/Decontam_process")
setwd("/Users/marlydmejia/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Desktop\ Documents-\ Patras\ Lab/Proposals/Manuscripts/Mine/2024_Rheumatoid_Arthritis_Vaginal/VSP80_2024updated/Decontam_process")
"
#IF NEEDED FOR MORE UPDATED ANALYSES
dir.create("decontam_items_final")
setwd("/Users/marlydmejia/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Desktop\ Documents-\ Patras\ Lab/VSP-RAproject_20240220/VSP80_2024updated/Decontam_process/decontam_items_final")
"
OTUmatrixfilt = phyloseq(otu_table(tFeatureTable,taxa_are_rows = FALSE), sample_data(metadata_set))

## ----see-meta-table--------------------------------------------------------
head(sample_data(OTUmatrixfilt))
quartz()

dffilt <- as.data.frame(sample_data(OTUmatrixfilt)) # Put sample_data into a ggplot-friendly data.frame
dffilt$LibrarySize <- sample_sums(OTUmatrixfilt)
dffilt <- dffilt[order(dffilt$LibrarySize),]
dffilt$Index <- seq(nrow(dffilt))
ggplot(data=dffilt, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(OTUmatrixfilt)$is.neg <- sample_data(OTUmatrixfilt)$Sample_or_Control == "Control"
# Make phyloseq object of presence-absence in negative controls and true samples
OTUmatrixfilt.pa <- transform_sample_counts(OTUmatrixfilt, function(abund) 1*(abund>0))
OTUmatrixfilt.pa.neg <- prune_samples(sample_data(OTUmatrixfilt.pa)$Sample_or_Control == "Control", OTUmatrixfilt.pa)
OTUmatrixfilt.pa.pos <- prune_samples(sample_data(OTUmatrixfilt.pa)$Sample_or_Control == "True_Sample", OTUmatrixfilt.pa)


## ----prevalence (0.1?)------------------------------------------------------------
contamdffilt.prev <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg")
table(contamdffilt.prev$contaminant)
head(which(contamdffilt.prev$contaminant))
# Make data.frame of prevalence in positive and negative samples
dffilt.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                        contaminant=contamdffilt.prev$contaminant)
ggplot(data=dffilt.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

contamdffilt.prev01b <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.1, batch="Batch", batch.combine = "minimum")
table(contamdffilt.prev01b$contaminant)
dffilt01b.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg), contaminant=contamdffilt.prev01b$contaminant)
ggplot(data=dffilt01b.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.1b)") + ylab("Prevalence (True Samples)")


#removed "batch= "Batch", batch.combine = "minimum" since read in as some batches having zero or one sample

## ----see-prev-.2-----------------------------------------------------------
contamdffilt.prev02 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.2, batch="Batch", batch.combine = "minimum")
table(contamdffilt.prev02$contaminant)
dffilt02.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg), contaminant=contamdffilt.prev02$contaminant)
ggplot(data=dffilt02.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.2)") + ylab("Prevalence (True Samples)")

## ----see-prev-.3-----------------------------------------------------------
contamdffilt.prev03 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.3, batch="Batch", batch.combine = "minimum")
table(contamdffilt.prev03$contaminant)
dffilt03.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg), contaminant=contamdffilt.prev03$contaminant)
ggplot(data=dffilt03.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.3)") + ylab("Prevalence (True Samples)")

## ----see-prev-.4-----------------------------------------------------------
contamdffilt.prev04 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.4, batch="Batch", batch.combine = "minimum")
table(contamdffilt.prev04$contaminant)
dffilt04.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg), contaminant=contamdffilt.prev04$contaminant)
ggplot(data=dffilt04.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.4)") + ylab("Prevalence (True Samples)")

## ----see-prev-.45-----------------------------------------------------------
contamdffilt.prev04_5 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.45, batch="Batch", batch.combine = "minimum")
table(contamdffilt.prev04_5$contaminant)
dffilt04_5.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg), contaminant=contamdffilt.prev04_5$contaminant)
ggplot(data=dffilt04_5.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.45)") + ylab("Prevalence (True Samples)")

## ----prevalence-.5---------------------------------------------------------
contamdffilt.prev05 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.5, batch="Batch", batch.combine = "minimum")
table(contamdffilt.prev05$contaminant)
dffilt05.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg), contaminant=contamdffilt.prev05$contaminant)
ggplot(data=dffilt05.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.5)") + ylab("Prevalence (True Samples)")


## ----prevalence-.05---------------------------------------------------------
contamdffilt.prev005 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.05, batch="Batch", batch.combine = "minimum")
table(contamdffilt.prev005$contaminant)
dffilt005.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg), contaminant=contamdffilt.prev005$contaminant)
ggplot(data=dffilt005.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.05)") + ylab("Prevalence (True Samples)")

#make sure you are in the right directory
write.csv(contamdffilt.prev01b,"contamdffilt.prev01b.csv", row.names = TRUE)
write.csv(contamdffilt.prev02,"contamdffilt.prev02.csv", row.names = TRUE)
write.csv(contamdffilt.prev03,"contamdffilt.prev03.csv", row.names = TRUE)
write.csv(contamdffilt.prev05,"contamdffilt.prev05.csv", row.names = TRUE)
write.csv(contamdffilt.prev005,"contamdffilt.prev005.csv", row.names = TRUE)

write.csv(dffilt01b.pa,"dffilt01b.pa.csv", row.names = TRUE)
write.csv(dffilt02.pa,"dffilt02.pa.csv", row.names = TRUE)
write.csv(dffilt03.pa,"dffilt03.pa.csv", row.names = TRUE)
write.csv(dffilt05.pa,"dffilt05.pa.csv", row.names = TRUE)
write.csv(dffilt005.pa,"dffilt005.pa.csv", row.names = TRUE)


#combine above tables to match everything by Feature ID
taxa_dffilt05 <- merge(dffilt05.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt05)[1] <- "Feature_ID"
taxa_contam05 <- merge(contamdffilt.prev05, taxonomy, by= 'row.names')
colnames(taxa_contam05)[1] <- "Feature_ID"
taxa_contam_dffilt05 <- merge(taxa_dffilt05, taxa_contam05, by= 'Feature_ID')
write.csv(taxa_contam_dffilt05,"taxa_contam_dffilt05.csv", row.names = FALSE)

taxa_dffilt04 <- merge(dffilt04.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt04)[1] <- "Feature_ID"
taxa_contam04 <- merge(contamdffilt.prev04, taxonomy, by= 'row.names')
colnames(taxa_contam04)[1] <- "Feature_ID"
taxa_contam_dffilt04 <- merge(taxa_dffilt04, taxa_contam04, by= 'Feature_ID')
write.csv(taxa_contam_dffilt04,"taxa_contam_dffilt04.csv", row.names = FALSE)

taxa_dffilt03 <- merge(dffilt03.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt03)[1] <- "Feature_ID"
taxa_contam03 <- merge(contamdffilt.prev03, taxonomy, by= 'row.names')
colnames(taxa_contam03)[1] <- "Feature_ID"
taxa_contam_dffilt03 <- merge(taxa_dffilt03, taxa_contam03, by= 'Feature_ID')
write.csv(taxa_contam_dffilt03,"taxa_contam_dffilt03.csv", row.names = FALSE)

taxa_dffilt02 <- merge(dffilt02.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt02)[1] <- "Feature_ID"
taxa_contam02 <- merge(contamdffilt.prev02, taxonomy, by= 'row.names')
colnames(taxa_contam02)[1] <- "Feature_ID"
taxa_contam_dffilt02 <- merge(taxa_dffilt02, taxa_contam02, by= 'Feature_ID')
write.csv(taxa_contam_dffilt02,"taxa_contam_dffilt02.csv", row.names = FALSE)

taxa_dffilt01b <- merge(dffilt01b.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt01b)[1] <- "Feature_ID"
taxa_contam01b <- merge(contamdffilt.prev01b, taxonomy, by= 'row.names')
colnames(taxa_contam01b)[1] <- "Feature_ID"
taxa_contam_dffilt01b <- merge(taxa_dffilt01b, taxa_contam01b, by= 'Feature_ID')
write.csv(taxa_contam_dffilt01b,"taxa_contam_dffilt01b.csv", row.names = FALSE)

taxa_dffilt005 <- merge(dffilt005.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt005)[1] <- "Feature_ID"
taxa_contam005 <- merge(contamdffilt.prev005, taxonomy, by= 'row.names')
colnames(taxa_contam005)[1] <- "Feature_ID"
taxa_contam_dffilt005 <- merge(taxa_dffilt005, taxa_contam005, by= 'Feature_ID')
write.csv(taxa_contam_dffilt005,"taxa_contam_dffilt005.csv", row.names = FALSE)

#manually removing columns of decontam-flagged taxa minus the 5 (code removed wrong taxa. uploading other csv returned incorrect variables - all characters)
#also, taxonomy uploaded where '-' became '.' so not understood is uploaded back to qiime2 
#------------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%% Qiime2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ########## FILTER OUT SPECIFIC TAXA ###########
#cd /Users/marlydmejia/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Desktop\ Documents-\ Patras\ Lab/Proposals/Manuscripts/Mine/2024_Rheumatoid_Arthritis_Vaginal

#Convert to qza file
biom convert -i VSP80_2024updated/exported/FULL_VSP80_gg2table/OTU_table_decontam.txt -o VSP80_2024updated/exported/FULL_VSP80_gg2table/OTU_feattable_decontam.biom --table-type="OTU table" --to-json

mkdir VSP80_2024updated/Decontam_process/decontam_items_final

#re-import new table to qiime2 for filtering of known contaminants
qiime tools import \
--input-path VSP80_2024updated/exported/FULL_VSP80_gg2table/OTU_feattable_decontam_JSN.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path VSP80_2024updated/Decontam_process/decontam_items_final/OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table VSP80_2024updated/Decontam_process/decontam_items_final/OTU_filtdecontam.qza \
--i-taxonomy VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales_650611; f__Pseudomonadaceae; g__Pseudomonas_E_647464; s__Pseudomonas_E_647464 veronii" \
--o-filtered-table VSP80_2024updated/VSP80_filtfeattable_gg2.qza

###################### Working files
qiime tools export \
--input-path VSP80_2024updated/VSP80_filtfeattable_gg2.qza \
--output-path VSP80_2024updated/exported/VSP80_filtfeattable_gg2
biom convert -i VSP80_2024updated/exported/VSP80_filtfeattable_gg2/feature-table.biom -o VSP80_2024updated/exported/VSP80_filtfeattable_gg2/feature-table.txt --to-tsv

qiime feature-table relative-frequency \
--i-table VSP80_2024updated/VSP80_filtfeattable_gg2.qza \
--o-relative-frequency-table VSP80_2024updated/VSP80_filtfeattable_gg2_relfreq
qiime tools export \
--input-path VSP80_2024updated/VSP80_filtfeattable_gg2_relfreq.qza \
--output-path VSP80_2024updated/exported/VSP80_filtfeattable_gg2_relfreq
biom convert -i VSP80_2024updated/exported/VSP80_filtfeattable_gg2_relfreq/feature-table.biom -o VSP80_2024updated/exported/VSP80_filtfeattable_gg2_relfreq/feature-table.txt --to-tsv



qiime feature-table filter-samples \
--i-table VSP80_2024updated/VSP80_filtfeattable_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--o-filtered-table VSP80_2024updated/Casecont2_feattable_gg2.qza
#removed VSP007 for abx use within week (kept VSP055 although has different circumstances - can always go back to file that was imported from decontam)

qiime feature-table summarize \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--o-visualization VSP80_2024updated/visualization/Casecont2_gg2.qzv \
--m-sample-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt

qiime tools export \
--input-path VSP80_2024updated/Casecont2_feattable_gg2.qza \
--output-path VSP80_2024updated/exported/Casecont2_feattable_gg2

biom convert -i VSP80_2024updated/exported/Casecont2_feattable_gg2/feature-table.biom -o VSP80_2024updated/exported/Casecont2_feattable_gg2/feature-table.txt --to-tsv


qiime feature-table relative-frequency \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--o-relative-frequency-table VSP80_2024updated/Casecont2_relfreq

qiime feature-table summarize \
--i-table VSP80_2024updated/Casecont2_relfreq.qza \
--o-visualization VSP80_2024updated/visualization/Casecont2_relfreq.qzv \
--m-sample-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt

qiime tools export \
--input-path VSP80_2024updated/Casecont2_relfreq.qza \
--output-path VSP80_2024updated/exported/Casecont2_feattable_gg2_relfreq

biom convert -i VSP80_2024updated/exported/Casecont2_feattable_gg2_relfreq/feature-table.biom -o VSP80_2024updated/exported/Casecont2_feattable_gg2_relfreq/feature-table.txt --to-tsv

              ####


qiime taxa barplot \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--i-taxonomy VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--o-visualization VSP80_2024updated/Visualization/barplot_Casecont2.qzv

qiime taxa collapse \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--i-taxonomy VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza \
--p-level 7 \
--o-collapsed-table VSP80_2024updated/Casecont2_feattable_species_gg2.qza

qiime composition add-pseudocount \
--i-table VSP80_2024updated/Casecont2_feattable_species_gg2.qza \
--o-composition-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza

qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column RA_status \
--o-visualization VSP80_2024updated/Visualization/compRA_Casecont2_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column Clinic \
--o-visualization VSP80_2024updated/Visualization/compClinic_Casecont2_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column RA_menopausebin \
--o-visualization VSP80_2024updated/Visualization/compRAmenopausebin_Casecont2_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column Menopausal_bin \
--o-visualization VSP80_2024updated/Visualization/compmenopausebin_Casecont2_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column ACPA_detectable \
--o-visualization VSP80_2024updated/Visualization/compACPAmerge_Casecont2_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column RF_detectable \
--o-visualization VSP80_2024updated/Visualization/compRFmerge_Casecont2_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column C_reactive_detectable \
--o-visualization VSP80_2024updated/Visualization/compCreacmerge_Casecont2_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column ACPA_disease \
--o-visualization VSP80_2024updated/Visualization/compACPAsplit_Casecont2_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column RF_disease \
--o-visualization VSP80_2024updated/Visualization/compRFsplit_Casecont2_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column C_reactive_disease \
--o-visualization VSP80_2024updated/Visualization/compCreacsplit_Casecont2_feattable_species_gg2.qzv

qiime feature-table filter-samples \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-where "[RA_Status]='Rheumatoid Arthritis'" \
--o-filtered-table VSP80_2024updated/CasesRA1_feattable_gg2.qza

qiime taxa collapse \
--i-table VSP80_2024updated/CasesRA1_feattable_gg2.qza \
--i-taxonomy VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza \
--p-level 7 \
--o-collapsed-table VSP80_2024updated/CasesRA_feattable_species_gg2.qza

qiime composition add-pseudocount \
--i-table VSP80_2024updated/CasesRA_feattable_species_gg2.qza \
--o-composition-table VSP80_2024updated/comp_CasesRA_feattable_species_gg2.qza

qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column Disease_Activity \
--o-visualization VSP80_2024updated/Visualization/compDiseaseAct_CasesRA_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column Clinic \
--o-visualization VSP80_2024updated/Visualization/compClinic_CasesRA_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column On_medication \
--o-visualization VSP80_2024updated/Visualization/compMeds_CasesRA_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column ACPA_detectable \
--o-visualization VSP80_2024updated/Visualization/compACPAdetect_CasesRA_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column RF_detectable \
--o-visualization VSP80_2024updated/Visualization/compRFdetect_CasesRA_feattable_species_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_species_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column C_reactive_detectable \
--o-visualization VSP80_2024updated/Visualization/compCreactdetect_CasesRA_feattable_species_gg2.qzv

STOPPED 032824 and 4/8/24 and 8/12/24 and 8/15/24 and 11/29/24









#genus level

qiime taxa collapse \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--i-taxonomy VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza \
--p-level 6 \
--o-collapsed-table VSP80_2024updated/Casecont2_feattable_genus_gg2.qza

qiime composition add-pseudocount \
--i-table VSP80_2024updated/Casecont2_feattable_genus_gg2.qza \
--o-composition-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza

qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column RA_status \
--o-visualization VSP80_2024updated/Visualization/compRA_Casecont2_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column Clinic \
--o-visualization VSP80_2024updated/Visualization/compClinic_Casecont2_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column Clinic_RA \
--o-visualization VSP80_2024updated/Visualization/compClinicRA_Casecont2_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column RA_menopausebin \
--o-visualization VSP80_2024updated/Visualization/compRAmenopausebin_Casecont2_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column Menopausal_bin \
--o-visualization VSP80_2024updated/Visualization/compmenopausebin_Casecont2_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column ACPA_detectable \
--o-visualization VSP80_2024updated/Visualization/compACPAmerge_Casecont2_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column RF_detectable \
--o-visualization VSP80_2024updated/Visualization/compRFmerge_Casecont2_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column C_reactive_detectable \
--o-visualization VSP80_2024updated/Visualization/compCreacmerge_Casecont2_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column ACPA_disease \
--o-visualization VSP80_2024updated/Visualization/compACPAsplit_Casecont2_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column RF_disease \
--o-visualization VSP80_2024updated/Visualization/compRFsplit_Casecont2_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_Casecont2_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column C_reactive_disease \
--o-visualization VSP80_2024updated/Visualization/compCreacsplit_Casecont2_feattable_genus_gg2.qzv

qiime feature-table filter-samples \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-where "[RA_Status]='Rheumatoid Arthritis'" \
--o-filtered-table VSP80_2024updated/CasesRA1_feattable_gg2.qza

qiime taxa collapse \
--i-table VSP80_2024updated/CasesRA1_feattable_gg2.qza \
--i-taxonomy VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza \
--p-level 6 \
--o-collapsed-table VSP80_2024updated/CasesRA_feattable_genus_gg2.qza

qiime composition add-pseudocount \
--i-table VSP80_2024updated/CasesRA_feattable_genus_gg2.qza \
--o-composition-table VSP80_2024updated/comp_CasesRA_feattable_genus_gg2.qza

qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column Disease_Activity \
--o-visualization VSP80_2024updated/Visualization/compDiseaseAct_CasesRA_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column Clinic \
--o-visualization VSP80_2024updated/Visualization/compClinic_CasesRA_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column On_medication \
--o-visualization VSP80_2024updated/Visualization/compMeds_CasesRA_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column ACPA_detectable \
--o-visualization VSP80_2024updated/Visualization/compACPAdetect_CasesRA_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column RF_detectable \
--o-visualization VSP80_2024updated/Visualization/compRFdetect_CasesRA_feattable_genus_gg2.qzv
qiime composition ancom \
--i-table VSP80_2024updated/comp_CasesRA_feattable_genus_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column C_reactive_detectable \
--o-visualization VSP80_2024updated/Visualization/compCreactdetect_CasesRA_feattable_genus_gg2.qzv

wget http://ftp.microbio.me/greengenes_release/current/2022.10.phylogeny.asv.nwk.qza
wget http://ftp.microbio.me/greengenes_release/current/VSP80_2024updated/2022.10.phylogeny.md5.nwk.qza
#Saved Phylogeny[Rooted] and didn't have to generate it myself

qiime tools export \
--input-path VSP80_2024updated/2022.10.phylogeny.md5.nwk.qza \
--output-path VSP80_2024updated/exported/rooted-tree_mergedstudy
#Exported 2022.10.phylogeny.md5.md5.qza as NewickDirectoryFormat to directory exported/rooted-tree_mergestudy


"
qiime diversity alpha-rarefaction \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--i-phylogeny VSP80_2024updated/2022.10.phylogeny.md5.nwk.qza \
--p-max-depth 10000 \
--p-metrics 'observed_features' \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-min-depth 7000 \
--p-steps 15 \
--p-iterations 50 \
--o-visualization VSP80_2024updated/Visualization/JK_otus_broad_allseq.qzv

qiime diversity alpha-rarefaction \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--i-phylogeny VSP80_2024updated/2022.10.phylogeny.md5.nwk.qza \
--p-max-depth 10000 \
--p-metrics 'shannon' \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-min-depth 7000 \
--p-steps 15 \
--p-iterations 50 \
--o-visualization VSP80_2024updated/Visualization/JK_shannon_broad_allseq.qzv
"

#STARTED 4/8/24 and 8/12/24 and 8/15/24
#make matrices: bray curtis
qiime diversity beta \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix VSP80_2024updated/Output/VSPMVScc1_BCDmatrix

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_BCDmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'RA_status' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RA_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_BCDmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'RA_status' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RA_bcdBeta_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_BCDmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'RA_menopausebin' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAmenopausebin_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_BCDmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'RA_menopausebin' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAmenopausebin_bcdBeta_permdisp


qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_BCDmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'Clinic_RA' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAclinic_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_BCDmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'Clinic_RA' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAclinic_bcdBeta_permdisp


#weighted normalized unifrac
qiime diversity beta-phylogenetic \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--i-phylogeny VSP80_2024updated/2022.10.phylogeny.md5.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-variance-adjusted FALSE \
--p-bypass-tips FALSE \
--o-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'RA_status' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RA_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'RA_status' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RA_wnuBeta_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'RA_menopausebin' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAmenopausebin_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'RA_menopausebin' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAmenopausebin_wnuBeta_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'Clinic_RA' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAclinic_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'Clinic_RA' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAclinic_wnuBeta_permdisp


############################added seqs to metadata files and samples-to-keep where needed
qiime diversity filter-distance-matrix \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_BCDmatrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--o-filtered-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--m-metadata-column 'Disease_Activity' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/CDAI_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--m-metadata-column 'Disease_Activity' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/CDAI_bcdBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--m-metadata-column 'ACPA_detectable' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/ACPA_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--m-metadata-column 'ACPA_detectable' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/ACPA_bcdBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--m-metadata-column 'Menopausal_Status' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/Menopause_gen_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--m-metadata-column 'Menopausal_Status' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/Menopause_gen_bcdBeta_permdisp

qiime diversity filter-distance-matrix \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--o-filtered-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--m-metadata-column 'Disease_Activity' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/CDAI_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--m-metadata-column 'Disease_Activity' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/CDAI_wnuBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--m-metadata-column 'ACPA_detectable' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/ACPA_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file Metadata_RAstudy_only.txt \
--m-metadata-column 'ACPA_detectable' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/ACPA_wnuBeta_permdisp





qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'Joint_space_narrowing' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAcont_jsn_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'Joint_space_narrowing' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAcont_jsn_wnuBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'Radiographic_Erosions' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAcont_re_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'Radiographic_Erosions' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAcont_re_wnuBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'AntiCCP_Ab_positive' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAcont_accp_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'AntiCCP_Ab_positive' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAcont_accp_wnuBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'RF_positive' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAcont_rf_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column 'RF_positive' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAcont_rf_wnuBeta_permdisp






qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'Joint_space_narrowing' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAjsn_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'Joint_space_narrowing' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAjsn_bcdBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'Joint_space_narrowing' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAjsn_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'Joint_space_narrowing' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAjsn_wnuBeta_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'jsn_aCCP' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAjsnCCP_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'jsn_aCCP' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAjsnCCP_bcdBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'jsn_aCCP' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAjsnCCP_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'jsn_aCCP' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAjsnCCP_wnuBeta_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'Radiographic_Erosions' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAre_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'Radiographic_Erosions' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAre_bcdBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'Radiographic_Erosions' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAre_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'Radiographic_Erosions' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAre_wnuBeta_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 're_aCCP' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAreCCP_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 're_aCCP' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAreCCP_bcdBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 're_aCCP' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAreCCP_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 're_aCCP' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAreCCP_wnuBeta_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'AntiCCP_Ab_positive' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAaccp_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'AntiCCP_Ab_positive' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAaccp_bcdBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'AntiCCP_Ab_positive' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAaccp_WNU_Beta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'AntiCCP_Ab_positive' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAaccp_WNU_Beta_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'RF_positive' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RARF_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'RF_positive' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RARF_bcdBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'RF_positive' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RARF_WNU_Beta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'RF_positive' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RARF_WNU_Beta_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'aCCP_RF' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAaCCPRF_bcdBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-BCdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'aCCP_RF' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAaCCPRF_bcdBeta_permdisp
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'aCCP_RF' \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAaCCPRF_wnuBeta_permanova
qiime diversity beta-group-significance \
--i-distance-matrix VSP80_2024updated/Output/RA-filtered-WNUdistance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_RA.txt \
--m-metadata-column 'aCCP_RF' \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization VSP80_2024updated/Visualization/RAaCCPRF_wnuBeta_permdis




%%%%alpha%%%%
  
DONE 04062024/08122024
qiime diversity alpha \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--p-metric observed_features \
--o-alpha-diversity VSP80_2024updated/Output/VSP80_alpha_otus
"qiime diversity alpha-correlation \
--i-alpha-diversity VSP80_2024updated/Output/VSP80_alpha_otus.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-method spearman \
--o-visualization VSP80_2024updated/Visualization/VSP80_alpha_otus_corr"
%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity VSP80_2024updated/Output/VSP80_alpha_otus.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--o-visualization VSP80_2024updated/Visualization/VSP80_alpha_otus_usethis

qiime diversity alpha \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--p-metric shannon \
--o-alpha-diversity VSP80_2024updated/Output/VSP80_alpha_shannon
"qiime diversity alpha-correlation \
--i-alpha-diversity VSP80_2024updated/Output/VSP80_alpha_shannon.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-method spearman \
--o-visualization VSP80_2024updated/Visualization/VSP80_alpha_shannon_corr"
%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity VSP80_2024updated/Output/VSP80_alpha_shannon.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--o-visualization VSP80_2024updated/Visualization/VSP80alpha_shannon_usethis
DONE 04062024/08122024/11242024








"
qiime diversity beta-rarefaction \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/VSP80_2024updated/2022.10.phylogeny.md5.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme RdGy_r \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_w_norm_unif_rarefaction_FULL.qzv \

qiime diversity beta-rarefaction \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--p-metric braycurtis \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme RdGy_r \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_bc_rarefaction_FULL.qzv \

qiime feature-table summarize \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--m-sample-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--o-visualization Visualization/Summary_Casecont2_featuretable_filt.qzv
#Saved Visualization to: Visualization/Summary_featuretableMM.qzv

%%bash -e
qiime diversity alpha-rarefaction \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--p-max-depth 1000 \
--p-metrics 'observed_features' \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-min-depth 5 \
--p-steps 25 \
--p-iterations 50 \
--o-visualization Visualization/JK_otu_max_longitudinal.qzv \

qiime diversity alpha-rarefaction \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--p-max-depth 1000 \
--p-metrics 'shannon' \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-min-depth 5 \
--p-steps 25 \
--p-iterations 50 \
--o-visualization Visualization/JK_shannon_max_longitudinal.qzv \
"





%%%%%%%%  - other diversity analyses - %%%%%%%%
#need to change file names accordingly


%%bash -e
qiime diversity pcoa \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_BCDmatrix.qza \
--p-number-of-dimensions 3 \
--o-pcoa VSP80_2024updated/Output/Casecont2_PCoA_BCD
#Saved PCoAResults to: Output/_PCoA_BCD.qza

%%bash -e
qiime emperor plot \
--i-pcoa VSP80_2024updated/Output/Casecont2_PCoA_BCD.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-ignore-missing-samples False \
--o-visualization VSP80_2024updated/Visualization/Casecont2_PCoA_BCD
#Saved Visualization to: Visualization/_PCoA_BCD.qzv

%%bash -e
qiime diversity pcoa-biplot \
--i-pcoa VSP80_2024updated/Output/Casecont2_PCoA_BCD.qza \
--i-features VSP80_2024updated/Casecont2_relfreq.qza \
--o-biplot VSP80_2024updated/Output/Casecont2_pcoa_BCD_biplot 
#Saved PCoAResults % Properties('biplot') to: Output/pcoa_BCD_biplot.qza

%%bash -e
qiime emperor biplot \
--i-biplot VSP80_2024updated/Output/Casecont2_pcoa_BCD_biplot.qza \
--m-sample-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-feature-metadata-file VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza \
--o-visualization VSP80_2024updated/Visualization/Casecont2_pcoaBCD_biplot
#Saved Visualization to: Visualization/pcoaBCD_biplot.qzv




%%bash -e
qiime diversity pcoa \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_WNUmatrix.qza \
--p-number-of-dimensions 3 \
--o-pcoa VSP80_2024updated/Output/Casecont2_PCoA_WNU
#Saved PCoAResults to: Output/_PCoA_WNUD.qza

%%bash -e
qiime emperor plot \
--i-pcoa VSP80_2024updated/Output/Casecont2_PCoA_WNU.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-ignore-missing-samples False \
--o-visualization VSP80_2024updated/Visualization/Casecont2_PCoA_WNU
#Saved Visualization to: Visualization/_PCoA_WNUD.qzv

%%bash -e
qiime diversity pcoa-biplot \
--i-pcoa VSP80_2024updated/Output/Casecont2_PCoA_WNU.qza \
--i-features VSP80_2024updated/Casecont2_relfreq.qza \
--o-biplot VSP80_2024updated/Output/Casecont2_pcoa_WNU_biplot 
#Saved PCoAResults % Properties('biplot') to: Output/pcoa_WNUD_biplot.qza

%%bash -e
qiime emperor biplot \
--i-biplot VSP80_2024updated/Output/Casecont2_pcoa_WNU_biplot.qza \
--m-sample-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-feature-metadata-file VSP80_2024updated/VSP80_01_gg2.tabletaxonomymd5.qza \
--o-visualization VSP80_2024updated/Visualization/Casecont2_pcoaWNU_biplot
#Saved Visualization to: Visualization/pcoaWNUD_biplot.qzv
Completed 081424 and 08/15/24









%%%%longitudinal%%%%  
  %%bash -e
qiime diversity filter-distance-matrix \
--i-distance-matrix VSP80_2024updated/Output/VSPMVScc1_BCDmatrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-where "mouse_line='HMB'" \
--o-filtered-distance-matrix Output/Casecont2_filtered-distance-matrix.qza
#Saved DistanceMatrix to: filtered-distance-matrix.qza

%%bash -e
qiime longitudinal first-distances \
--i-distance-matrix Output/Casecont2_filtered-distance-matrix.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-state-column tp_days \
--p-individual-id-column mouse_origid \
--p-replicate-handling error \
--o-first-distances Output/Casecont2_BC_FirstDistances 
#Saved SampleData[FirstDifferences] to: _BC_FirstDistances.qza

qiime tools export \
--input-path Output/Casecont2_BC_FirstDistances.qza \
--output-path exported/Casecont2_BC_FirstDistances

qiime tools export \
--input-path Output/Casecont2_pcoa_BCD_biplot.qza \
--output-path exported/Casecont2_pcoa_BCD_biplot.qza

qiime tools export \
--input-path VSP80_2024updated/Output/VSPMVScc1_BCDmatrix.qza \
--output-path exported/Casecont2_BCDmatrix.qza


qiime diversity pcoa \
--i-distance-matrix Output/Unrarefied_wnu_matrices_longitudinal.qza \
--p-number-of-dimensions 3 \
--o-pcoa Output/Casecont2_PCoA_wnu
qiime emperor plot \
--i-pcoa Output/Casecont2_PCoA_wnu.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-ignore-missing-samples False \
--o-visualization Visualization/Casecont2_PCoA_wnu
%%bash -e
qiime diversity pcoa-biplot \
--i-pcoa Output/Casecont2_PCoA_wnu.qza \
--i-features Output/Casecont2_relative.qza \
--o-biplot Output/Casecont2_pcoa_wnu_biplot
%%bash -e
qiime emperor biplot \
--i-biplot Output/Casecont2_pcoa_wnu_biplot.qza \
--m-sample-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-feature-metadata-file mergedstudy.gg2.tabletaxonomy.qza \
--o-visualization Visualization/Casecont2_pcoawnu_biplot

qiime feature-table filter-samples \
--i-table VSP80_2024updated/Casecont2_feattable_gg2.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--p-where "[body_site]='vaginal tract'" \
--o-filtered-table HMB_longitudinal-table.qza

qiime longitudinal first-distances \
--i-distance-matrix Output/Casecont2_filtered-distance-matrix.qza \
--m-metadata-file HMBpaper_metadata_longitudinalsubset_noNAN.txt \
--p-state-column tp_consecswabs \
--p-individual-id-column mouse_throughstudies \
--p-replicate-handling error \
--o-first-distances Output/Casecont2_BC_FirstDistances_staging 

qiime tools export \
--input-path Output/Casecont2_BC_FirstDistances_staging.qza \
--output-path exported/Casecont2_BC_FirstDistances_staging





qiime diversity beta-phylogenetic \
--i-phylogeny Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/VSP80_2024updated/2022.10.phylogeny.md5.nwk.qza \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/VSP80_2024updated/Casecont2_feattable_gg2.qza \
--p-metric weighted_normalized_unifrac \
--p-variance-adjusted TRUE \
--o-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/
  

"#Longitudinal"
qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_longitudinal.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column cohort_number \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_longCohortcombnum_wnupermanova
#Saved Visualization to: Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_longitudinal.qza \
--m-metadata-file VSP80_2024updated/Metadata_RAstudy_master_categorized_v2excess.txt \
--m-metadata-column cohort_number \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_longCohortcombnum_wnupermdisp

qiime tools export \
--input-path /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_longCohortcombnum_wnupermanova.qzv \
--output-path  /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/HMBfinal_longCohortcombnum_wnupermanova
qiime tools export \
--input-path /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_longCohortcombnum_wnupermdisp.qzv \
--output-path  /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/HMBfinal_longCohortcombnum_wnupermdisp





#------------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%% Valencia clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  ###################### Assigning CST!!!
cd VSP80_2024updated/Valencia\ copy/  
python3 /path/to/convert_qiime.py \
VSP_MVS_01_gg2.tabletaxonomymd5/taxonomy.csv \
OTU_gg2taxatable_Valencia.csv

cd VSP80_2024updated/Valencia
python3 Valencia.py -r CST_centroids_012920.csv -i OTU_gg2taxatable_Valencia.csv -o Valencia_classifications_20240813_VSP80 -p Valencia_classifications_plot_20240813_VSP80

"It is imperative that phylotype column headings match those used by the VALENCIA reference centroids which generally take the form of
“phylotype rank underscore phylotype name” (e.g., g_Bifidobacterium). All phylotypes should be summarized to the genus rank or higher except for the following:
Lactobacillus spp., Gardnerella spp., Prevotella spp., Atopobium spp., Sneathia spp, Mobiluncus spp.
These key phylotypes appear as “Genus underscore species” (e.g., Lactobacillus_crispatus, Gardnerella_vaginalis)"


"VALENCIA has the follow required arguments:

-ref, --reference : path to the reference centroids file (provided)

-i, --input : CSV file of the sample dataset with column 1 named “sampleID” containing unique sample names and column 2 named “read_count” containing each sample’s total read count. The remaining columns are expected to be taxa read counts, with the taxa name as the header.

-o, --output : File prefix to store the output, all of the added CST information will appear as the final columns of the dataset

And the following optional arguments:

-p, --plot : File prefix to store diagnostic plot in, not generated unless this argument is provided"

#python3 /path/to/Valencia.py -ref /path/to/CST_profiles_012920.csv -i /path/to/test_dataset.csv -o /path/to/test_out -p /path/to/test




