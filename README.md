# Splicing-Pipeline

**1. Mapping using STARsolo**
-
We used [STAR](https://github.com/alexdobin/STAR/tree/master) as the primary mapper for this pipeline. For this step you will need R1 and R2 fastq files as well as the GTF file and the index built using STAR. The raw fastq files have been procced with STARsolo using `--soloFeatures SJ` to get the counts for annotated and novel splice junctions. The output of shold be present in the `SJ.out.tab` file. The outputh should be like this:
```
chr1    3109059 3226649 1       1       0       0       4       43
chr1    3122081 3133046 1       1       0       1       0       15
chr1    3145384 3146420 2       2       0       3       0       42
chr1    3220998 3226168 2       2       0       0       1       22
chr1    3271768 3271922 0       0       0       3       0       45
chr1    3277541 3283661 2       2       1       6       0       40
chr1    3281209 3283661 2       2       0       5       0       43
chr1    3287192 3491924 2       2       1       2       0       23
```

The columns have the following meaning:
```
column 1: chromosome
column 2: first base of the intron (1-based)
column 3: last base of the intron (1-based)
column 4: strand (0: undefined, 1: +, 2: -)
column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
column 6: 0: unannotated, 1: annotated in the splice junctions database. Note that in 2-pass mode, junctions detected in the 1st pass are reported as annotated, in addition to annotated junctions from GTF.
column 7: number of uniquely mapping reads crossing the junction
column 8: number of multi-mapping reads crossing the junction
column 9: maximum spliced alignment overhang
```

We are useing the raw tab file for the junction abundace and will apply ad-hoc filtartion troughout the pipeline. For more information about how to build genome index and mapping using STAR please take a look at examples codes.


Here's a simple flowchart

```mermaid
graph TD;
fastqs-->STARsolo
STARsolo-->Splicing-pipeline
Splicing-pipeline-->GEDI-PM

```


**2. Grouping the junction abundance**
--

We grouped splice junctions according to the coordinates of their intronic regions into “local junction variants” (LJVs). An LJV is defined as a set of junctions that share either the first or the last coordinate.

For each junction within an LJV, we extracted the junction abundance count for every cell, which we refer to as the **M1 count**. The **M2 count** for a given junction in each cell is defined as the sum of the counts of its alternative junctions—that is, all other junctions within the same LJV. In other words, for each junction and cell, the number of reads supporting that junction is contrasted with the total number of reads supporting its alternative junctions. For example, if $( M1_{ij} )$ represents the count for junction $j$ in cell $i$, then the M2 count can be written as:


```math
M2_{ij} = \sum_{\substack{k \in J \\ k \neq j}} M1_{ik}
```


where $J$ denotes the set of all junctions in the LJV. The grouping method is determined by whether the first or the last coordinate is used, with an appended `-E` or `-S` added to the event IDs accordingly. This approach may result in a single junction receiving two different measurements in the M2 counts (in the inclusion matrix) while having identical measurements in the M1 matrix.

We also handle sample-specific junctions. If a junction is present in only a subset of samples, a corresponding vector of zeros is applied to the M1 matrix for the samples in which the junction is absent, and the M2 measurements are computed accordingly. This figure shows a perfect examples of two LJVs.

![alt text](https://github.com/Arshammik/Splicing-Pipeline/blob/main/Markdown_figure.jpg?raw=true)


Subsequently, we employed a Poisson read splitting (thinning) technique to divide the raw counts into train and test data sets to address the double-dipping problem (wherein the same data set is used to generate and test a hypothesis). This was followed by a filtration of highly variable features in both modalities. 


# Usage

The first step of the pipeline is to make a list of junction abundance and event data **per sample**. `eventdata` has reports for each junction and has the same number of rows as M1 and M2 matricies contaies the chromosome number, start and end coordinates as well as other quality and grouping metrics. This step could be done using `multigedi_make_junction_ab`.
### 1. `multigedi_make_junction_ab`

**Example:**

```r
SJ_object <- multigedi_make_junction_ab(STARsolo_SJ_dirs= c("./example_star_solo_outout/Solo.out/SJ/"),
                                        sample_ids=c("SMP_1"),
                                        use_internal_whitelist = TRUE)
```

**Arguments:**


- `STARsolo_SJ_dirs`   A character vector or list of strings representing the paths to STARsolo SJ directories. Each directory should contain the raw splicing junction output files.
sample_ids

- `sample_ids`   A character vector or list of unique sample IDs corresponding to each directory in `STARsolo_SJ_dirs`.(It will attached to the barcodes in the final matrices).

- `use_internal_whitelist`   A logical flag (default `TRUE`) indicating whether to use the internal STARsolo whitelist located at `../Gene/filtered/barcodes.tsv` for each sample when `white_barcode_lists` is `NULL`.

- `white_barcode_lists`   A list of character vectors, each containing barcode whitelist(s) for the corresponding sample. If `NULL` (default), the function uses the internal STARsolo whitelist if `use_internal_whitelist` is `TRUE`.



After getting the object with sample specefic junction abundance matrices and related event data you need to use `multigedi_make_m1` to get all in one M1 and eventdata.

### 2. `multigedi_make_m1`

**Example:**

```r
m1_obj <- multigedi_make_m1(junction_ab_object= SJ_object)

# # the out puth should be a list with two objects, one is M1 matrix and the other one should be eventdata
# summary(m1_obj)
# 
# >                     Length Class      Mode
# > m1_inclusion_matrix 663320 dgCMatrix  S4  
# > event_data              13 data.table list
```

**Arguments:**
- `junction_ab_object`   A named list that function `multigedi_make_junction_ab` returns where each element represents a sample's junction abundance data. Each element must contain eventdata` and a sparse matrix

After obtaining M1 and the event data, you must generate an exclusion matrix, M2. If you intend to perform a count split, ensure that this operation is done prior to creating M2. This order is crucial because it guarantees that the resulting `m1_test` and `m1_train` subsets remain statistically independent. In this context, you should use the `multigedi_countsplit` function, which is essentially an integrated version of the main function from the [countsplit R package](https://github.com/anna-neufeld/countsplit/blob/develop/R/countsplit.R).

### 3. `multigedi_countsplit`

**Example:**

```r
m1 <- m1_obj$m1_inclusion_matrix
multigedi_countsplit(m1_inclusion_matrix= m1,
                    folds = 2,
                    epsilon  = c(0.5, 0.5),
                    object_names = "m1")
```
**Arguments:**

- `m1_inclusion_matrix`   A dense or sparse  numerical matrix to be split.
- `folds`   An integer specifying how many folds you would like to split your data into.
- `epsilon`   epsilon A vector, which has length `folds`, that stores non-zero elements that sum to one. Determines the proportion of information from X that is allocated to each fold. Defult is `c(0.5, 0.5)`.
- `object_names`   A character string specifying the base name for output train/test objects. The deafult is "m1".

There is no need to assign a variable for this command; it will automatically store the output in the `.GlobalEnv` as `object_names_test` and `object_names_train`.<br/>
After performing the count splitting, derive the M2 matrix from the M1 matrix. This step should be executed only once for each individual M1 matrix, regardless of whether it is original, test, or training data.

### 4. `multigedi_make_m2`

**:warning:** This step may be memory intensive. Please ensure that an adequate amount of memory is allocated based on the size of the M1 matrix.


**Example:**

```r
m2_test <- multigedi_make_m2(m1_inclusion_matrix= m1_test, eventdata=m1_obj$event_data)
m2_train <- multigedi_make_m2(m1_inclusion_matrix= m1_train, eventdata=m1_obj$event_data)
```

**Arguments:**

- `m1_inclusion_matrix`   A sparse matrix to be modified and used for creating the M2 matrix.
- `eventdata`   A data table containing event information with at least `group_id` and an index column.


At this point, the pipeline for obtaining junction abundance measurements is complete, and you can calculate PSI values using the GEDI pan-modal tool from the M1 and M2 matrices. In addition, if you want to integrate gene expression with splicing data, it is recommended to generate the gene expression count matrix directly from the STARsolo outputs using the `multigedi_make_gene` function. This approach ensures consistency by relying on a unified mapping strategy.


### 5. `multigedi_make_gene`

**Example:**


```r
gene_expression <- multigedi_make_gene(expression_dirs = c("./example_star_solo_outout/Solo.out/Gene/"),
                                       sample_ids=c("SMP_1"),
                                       use_internal_whitelist=TRUE)
```
**Arguments:**

- `expression_dirs` A character vector or list of strings representing paths to directories containing the gene expression matrix (`matrix.mtx`), barcodes, and features files.

- `sample_ids`   A character vector or list of unique sample IDs corresponding to each directory in `expression_dirs`.(It will attached to the barcodes in the final matrices).

- `use_internal_whitelist`   A logical flag (default `TRUE`) indicating whether to use the internal STARsolo whitelist located at `../Gene/filtered/barcodes.tsv` for each sample when `white_barcode_lists` is `NULL`.

- `white_barcode_lists`   A list of character vectors, each containing barcode whitelist(s) for the corresponding sample. If `NULL` (default), the function uses the internal STARsolo whitelist if `use_internal_whitelist` is `TRUE`. Note that for the external baecodes you need to remove the number that denotes the GEM well, meaning that the `AAACCCAAGGAGAGTA-1` should turn into `AAACCCAAGGAGAGTA` for effective filtration.



If STARsolo mapping was run with the option `--soloFeatures Velocyto`, the resulting output should have a structure similar to the following:

```bash
.
├── SJ.out.tab
└── Solo.out
    ├── Barcodes.stats
    ├── Gene
    ├── GeneFull
    ├── SJ
    └── Velocyto
```

Using the function `multigedi_make_velo`, you can generate inclusion and exclusion matrices for intronic and exonic reads at the gene feature level. These matrices serve a similar purpose to the M1 and M2 matrices for junction abundance.

### 5. `multigedi_make_velo`
**Example:**

```r
velocyto_obj <- multigedi_make_velo(velocyto_dirs=c("./example_star_solo_outout/Solo.out/Velocyto/"), 
                    sample_ids=c("SMP_1"),
                    use_internal_whitelist = TRUE,
                    merge_counts = TRUE)

# > summary(velocyto_obj)
# >          Length    Class     Mode
# > spliced   322344680 dgCMatrix S4  
# > unspliced 322344680 dgCMatrix S4 
```
**Arguments:**

- `velocyto_dirs`   A character vector or list of strings representing paths to Velocyto directories. Each directory should contain subdirectories (`filtered` or `raw`) with required files.

- `sample_ids`   A character vector or list of unique sample IDs corresponding to each directory in `velocyto_dirs`.

- `whitelist_barcodes` A list of character vectors, each containing barcode whitelist(s) for the corresponding sample. If `NULL` (default), the function uses the filtered barcodes file if `use_internal_whitelist` is `TRUE`.

- `use_internal_whitelist` A logical flag (default `TRUE`) indicating whether to use the `filtered` data for barcode filtration. If `FALSE`, the `raw` data is used, and no barcode filtration is applied unless a whitelist is provided.
- `merge_counts`: A logical flag that controls whether spliced and unspliced count matrices from all samples should be combined.
  
  - **Default:** `FALSE`
  
  - **When set to `TRUE`:**
    - All spliced matrices are merged into a single large spliced matrix.
    - All unspliced matrices are merged into a single large unspliced matrix.
  
  - **When set to `FALSE`:**
    - The function returns a list where each sample is represented by two `dgCMatrix` objects: one for spliced counts and one for unspliced counts.


