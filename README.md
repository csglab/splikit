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


**2. Grouping the junction abundance**
--
# Revised Document on Splice Junction Analysis

We grouped splice junctions according to the coordinates of their intronic regions into “local junction variants” (LJVs). An LJV is defined as a set of junctions that share either the first or the last coordinate.

For each junction within an LJV, we extracted the junction abundance count for every cell, which we refer to as the **M1 count**. The **M2 count** for a given junction in each cell is defined as the sum of the counts of its alternative junctions—that is, all other junctions within the same LJV. In other words, for each junction and cell, the number of reads supporting that junction is contrasted with the total number of reads supporting its alternative junctions.

## Caution on Summation in LaTeX

When expressing the summation of matrices or counts in LaTeX, ensure that the indices and limits of summation are clearly defined to avoid ambiguity. For example, if \( M1_{ij} \) represents the count for junction \( j \) in cell \( i \), then the M2 count can be written as:
\[
M2_{ij} = \sum_{\substack{k \in J \\ k \neq j}} M1_{ik},
\]
where \( J \) denotes the set of all junctions in the LJV. Such clarity is essential to correctly communicate the intended computation.

The grouping method is determined by whether the first or the last coordinate is used, with an appended `-E` or `-S` added to the event IDs accordingly. This approach may result in a single junction receiving two different measurements in the M2 counts (in the inclusion matrix) while having identical measurements in the M1 matrix.

We also handle sample-specific junctions. If a junction is present in only a subset of samples, a corresponding vector of zeros is applied to the M1 matrix for the samples in which the junction is absent, and the M2 measurements are computed accordingly.

 
 
Subsequently, motivated by (Anna Neufeld et.al)[https://arxiv.org/abs/2207.00554] we employed a Poisson read splitting (thinning) technique to divide the raw counts into train and test data sets to address the double-dipping problem (wherein the same data set is used to generate and test a hypothesis). This was followed by a filtration of highly variable features in both modalities. 

