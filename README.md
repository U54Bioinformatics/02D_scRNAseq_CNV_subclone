# InferCNV

### input files

*fibroblasts counts is used as a reference here.*

- filtered counts (the "filtered_counts" in the demo script): 

  the counts of the high quality cells in the format of tab-delimited plain text file, looking like so:  

| Gene.Symbol | Cell1 | Cell2 | ...  |
| ----------- | ----- | ----- | ---- |
| A           | 1     | 2     | ...  |
| B           | 0     | 5     | ...  |
| C           | 2     | 6     | ...  |
| ...         | ...   | ...   | ...  |

- reference counts: a plain text file looking like so:  

| Gene.Symbol | Ref1 | Ref2 | ...  |
| ----------- | ---- | ---- | ---- |
| C           | 1    | 0    | ...  |
| A           | 0    | 0    | ...  |
| B           | 0    | 1    | ...  |
| ...         | ...  | ...  | ...  |

- cell annotations (the "singler_annot" in the demo script): 

  a plain text file with column headers:  

```r
Cell annot
Cell1 Epithelial
Cell2 Endotheial
... ...
```

- gene position file (the "gene_order_file" in the "CreateInfercnvObject" function): 

  The format is tab-delimited and has no column header, simply providing the gene name, chromosome, and gene span:

```R
A    Chr2	14363	27239
B    ChrX 761586 762902
C    Chr5 1152288 1167411
...    ...
```

### scripts

```R
# this script is reproducible in the R docker image on unicron

# load librarys
library(infercnv)
library(dplyr)
library(data.table)

set.seed(41)
# import raw data
## counts
### filtered counts
filtered_counts <- 
	fread('xxx/xxxfiltered_counts.txt', stringsAsFactors = F, check.names =F) %>% 
		as.data.frame
### fibroblast counts (reference)
fibro_counts <- 
	fread('xxx/80_scrna_cnv_normalization/Fibroblasts.counts.txt', stringsAsFactors = F,
        check.names =F) %>% 
		as.data.frame
## singler annot
singler_annot <- 
	read.delim('xxx/singler_filtered_annot.txt', stringsAsFactors = F, check.names = F)

# data cleaning
## change geneid colum name to genesymbol
colnames(fibro_counts)[1] = "Gene.Symbol"
## generate and write the annotation file
fibro_annot <- data.frame(Cell = colnames(fibro_counts)[-1], annot = "fibroblast")
write.table(rbind(fibro_annot, singler_annot), 
            "filtered_fibro-ref.txt", sep = "\t", quote = F, 
            col.names = F, row.names = F)
## generate counts matrix
filtered_fibro_counts_matrix <- 
inner_join(filtered_counts, fibro_counts, by = "Gene.Symbol") %>% 
	`rownames<-`(c(.$Gene.Symbol)) %>% 
	.[, !colnames(.) %in% "Gene.Symbol"] %>% 
	as.matrix
 
# do infercnv
# creat infercvnobj
infercnvobj <- CreateInfercnvObject(raw_counts_matrix=filtered_fibro_counts_matrix, 
                                    annotations_file="filtered_fibro-ref.txt", delim="\t", 
                                    gene_order_file="xxx/gencode_v19_gene_pos.txt",
                                    ref_group_names=c("fibroblast"),  #an alternative: NULL
                                    chr_exclude = NULL #an alternative: c("ChrM", "ChrX")
                                   )

# perform infercnv operations to reveal cnv signal
infercnvobj = infercnv::run(infercnvobj, 
cutoff=0.1,
out_dir="./out-fibro-ref-k3",
cluster_by_groups=F, 
k_obs_groups = 3,
denoise=T, HMM=T,  
num_threads = 23)
```

