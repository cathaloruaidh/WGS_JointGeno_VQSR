# WGS_JointGeno_VQSR

# Dependencies
The main sections of the pipeline (modules 1 to 10) use the following programs: 
- `gatk` v3.8
- `bcftools` v1.11
- `VCFtools` v0.1.15
- `R` v3.6.0
- `htslib` v1.4.1 (for `bgzip` and `tabix`)

The remaining sections of the pipeline (modules 11 to 15) additionally use the following programs: 
- `plink` v1.90
- `peddy` v0.4.7
- `pdflatex` v3.1415926-2.5-1.40.14 (TeX Live 2013)
- `vt` v0.5772-60f436c3

The following resource files are required for VQSR:
- Mills/1000 Genomes  gold standard indels
- dbSNP 
- 1000 Genomes Phase 3 SNVs
- Omni 2.5 SNVs
- HapMap 3.3 SNVs

# Help
```
-b, --build - hg19 or GRCh38 [GRCh38]
-c, --clean - Remove all input files once finished [false]
-d, --dbsnp <FILE> - dbSNP VCF file  [v150]
-e, --end <INT> - Tool to finish on [15]
-f, --file <FILE> - Name of source file
-g, --gatk <FILE> - GATK file [GATK v3.8]
-h, --help - Output this message
-m, --mem-min <INT> - Java minimum memory [6G]
-o, --out <FILE> - Output file prefix ["Variant_Filter"]
-r, --reference <FILE> - FASTA Reference file [GRCh38]
-s, --start <INT> - Tool to start on [1]
-t, --threads <INT> - Number of threads for multithreaded processes [4]
-v, --verbose <INT> - Set verbosity level [3]
-x, --mem-max <INT> - Java maximum memory [6G]
    --snpTS <FLOAT> - Tranche Sensitivity for SNPs [99.9]
    --indelTS <FLOAT> - Tranche Sensitivity for indels; [99.0]
    --fam <FILE> - FAM file ; []
 ```
