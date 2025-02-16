# Splicing-Pipeline Summary
This pipeline provides a streamlined approach for processing and analyzing splicing data from single-cell RNA sequencing. It begins with mapping raw reads using **STARsolo**, capturing splice junctions in the `SJ.out.tab` file. Next, junctions sharing intron coordinates are grouped into Local Junction Variants (LJVs), and counts are split into **M1** (junction-specific reads) and **M2** (reads supporting alternative junctions). A Poisson-based thinning method helps split these counts into independent train/test datasets, ensuring unbiased downstream analyses. Additional tools within the pipeline facilitate generating gene expression matrices, combining splicing and gene-level measurements, and preparing Velocyto (spliced/unspliced) data, offering a comprehensive workflow for single-cell splicing quantification.

![alt text](https://github.com/Arshammik/Splicing-Pipeline/blob/main/Markdown_figure.jpg?raw=true)


**1. Mapping with STARsolo**  
- Map paired-end fastq files using STARsolo with a GTF and pre-built index.  
- Generates `SJ.out.tab` containing junction info (chr, start, end, strand, motif, read counts).

**2. Junction Grouping & Count Splitting**  
- Group junctions into Local Junction Variants (LJVs) sharing an intron coordinate.  
- For each junction, define:
  - **M1**: Junction count.
  - **M2**: Sum of counts for alternative junctions:
    ```math
    M2_{ij} = \sum_{k \neq j} M1_{ik}
    ```
- Apply Poisson thinning to split counts into independent train/test datasets.

**3. Key Functions**  
- `multigedi_make_junction_ab`: Generate per-sample junction abundance.  
- `multigedi_make_m1`: Combine M1 matrices with event data.  
- `multigedi_countsplit`: Split M1 into train/test sets.  
- `multigedi_make_m2`: Create the M2 (exclusion) matrix.  
- `multigedi_make_gene` & `multigedi_make_velo`: Generate gene expression and Velocyto (spliced/unspliced) matrices.

**Visualization**

```mermaid
graph TD;
fastqs-->STARsolo;
STARsolo-->Splicing-pipeline;
Splicing-pipeline-->GEDI-PM;
