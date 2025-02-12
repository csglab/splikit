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

