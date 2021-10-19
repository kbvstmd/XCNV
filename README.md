# 1. Introduction
 X-CNV is a tool to predict CNV pathogenicity using an XGBoost classifier.<br><br>
 X-CNV calculates a meta-voting prediction (MVP) score to quantitatively evaluate disease-causing probability. It consists of the most comprehensive CNV data and annotations by integrating various publicly available genetic variant repositories. The features covering the genomics, genome region, variation types, and population genetics properties are taken into account to boost the prediction power. More importantly, a meta-voting prediction (MVP) score is proposed to measure the CNV pathogenic effect quantitatively, which can be used to determine the CNV pathogenicity.The reference genome version used by X-CNV is GRCh37/hg19. <br>
 Cite:Zhang L, Shi J, Ouyang J, Zhang R, Tao Y, Yuan D, Lv C, Wang R, Ning B, Roberts R, Tong W, Liu Z, Shi T. X-CNV: genome-wide prediction of the pathogenicity of copy number variations. Genome Med. 2021 Aug 18;13(1):132. doi: 10.1186/s13073-021-00945-4. PMID: 34407882.
# 2. Requirements
 The local version X-CNV requires two R packages, data.table and xgboost(your may get issues when you install the xgboost packages you can see the https://github.com/kbvstmd/XCNV/issues/1 for help), and Bedtools v2.26.0. If the R packages and bedtools cannot be installed automatically, users can install them manually. The executable file of bedtools should be placed in ./tools/. ## Memory limit >=8G
# 3. Installation
```bash
git clone https://github.com/kbvstmd/XCNV.git
cd XCNV
sh Install.sh
```
# 4. Usage and example
## Usage:
```bash
./bin/XCNV prefix.bed
```
The output filename: prefix.output.csv

## Example:
```bash
./bin/XCNV ./example_data/1.bed
```
The results can be seen in the 1.output.csv

# 5. Input & output
Input file format (The columns are separated by TAB key and the header is not required): <br><br>

2   2222999 3000222 gain <br>

Column 1: The chromosome (no “chr”) <br>
Column 2: Start <br>
Column 3: End <br>
Column 4: CNV type (gain or loss) <br>

The output file has 35 columns and is provided as Comma-Separated Values (CSV) format. <br>

<table>
    <thead>
    <th >Columns</th>
    <th >Description</th>
    <th>Category</th>
    </thead>
    <tr><td>Chr</td><td>Chromosome</td><td>Input</td></tr>
<tr><td>Start</td><td>Start position</td><td>Input</td></tr>
<tr><td>End</td><td>End position</td><td>Input</td></tr>
<tr><td>CNV type</td><td>CNV type (gain or loss)</td><td>Input</td></tr>
<tr><td>FATHMM score</td><td>FATHMM Dnase score for the CNV region</td><td>Coding</td></tr>
<tr><td>LR score</td><td>LR score for the CNV region</td><td>Coding</td></tr>
<tr><td>LRT score</td><td>LRT Dnase score for the CNV region</td><td>Coding</td></tr>
<tr><td>MutationAssessor score</td><td>MutationAssessor score for the CNV region</td><td>Coding</td></tr>
<tr><td>MutationTaster score</td><td>MutationTaster score for the CNV region</td><td>Coding</td></tr>
<tr><td>Polyphen2-HDIV score</td><td>Polyphen2_HDIV score for the CNV region</td><td>Coding</td></tr>
<tr><td>Polyphen2-HVAR score</td><td>Polyphen2_HVAR score for the CNV region</td><td>Coding</td></tr>
<tr><td>RadialSVM score</td><td>RadialSVM Dnase score for the CNV region</td><td>Coding</td></tr>
<tr><td>SIFT score</td><td>SIFT Dnase score for the CNV region</td><td>Coding</td></tr>
<tr><td>VEST3 score</td><td>VEST3 score for the CNV region</td><td>Coding</td></tr>
<tr><td>pLI</td><td>Probability of being loss-of-function intolerant</td><td>Coding</td></tr>
<tr><td>Episcore</td><td>A computational method to predict haploinsufficiency leveraging epigenomic data from a broad range of tissue and cell types by machine learning methods.</td><td>Coding</td></tr>
<tr><td>GHIS</td><td>An integrative approach to predicting haploinsufficient genes</td><td>Coding</td></tr>
<tr><td>CADD score</td><td>Average CADD score for the CNV region</td><td>Genome-wide</td></tr>
<tr><td>GERP</td><td>GERP++_RS Dnase score for the CNV region</td><td>Genome-wide</td></tr>
<tr><td>phyloP100way</td><td>phyloP100way_vertebrate score for the CNV region</td><td>Genome-wide</td></tr>
<tr><td>phyloP46way</td><td>phyloP46way_placental score for the CNV region</td><td>Genome-wide</td></tr>
<tr><td>SiPhy29way</td><td>SiPhy_29way_logOdds score for the CNV region</td><td>Genome-wide</td></tr>
<tr><td>cdts-1st</td><td>The coverage ratio between  CDTS percentile < 1% and the CNV region</td><td>Noncoding</td></tr>
<tr><td>cdts-5th</td><td>The coverage ratio between  CDTS percentile < 5% and the CNV region</td><td>Noncoding</td></tr>
<tr><td>pELS</td><td>The coverage of proximal enhancer-like sequence (pELS) within the CNV region</td><td>Noncoding</td></tr>
<tr><td>CTCF-bound</td><td>The coverage of CTCF-bound sequence within the CNV region</td><td>Noncoding</td></tr>
<tr><td>PLS</td><td>The coverage of promoter-like sequence within the CNV region</td><td>Noncoding</td></tr>
<tr><td>dELS</td><td>The coverage of distal enhancer-like sequence within the CNV region</td><td>Noncoding</td></tr>
<tr><td>CTCF-only</td><td>The coverage of CTCF-only sequence within the CNV region</td><td>Noncoding</td></tr>
<tr><td>DNase-H3K4me3</td><td>The coverage of DNase-H3K4me3 sequence within the CNV region</td><td>Noncoding</td></tr>
<tr><td>gain-PAF</td><td>Population allele frequency for duplication</td><td>Universal</td></tr>
<tr><td>Length</td><td>CNV length</td><td>Universal</td></tr>
<tr><td>loss-PAF</td><td>Population allele frequency for deletion</td><td>Universal</td></tr>
<tr><td>Type.1</td><td>CNV type (gain or loss code as 1 or 0)</td><td>Universal</td></tr>
<tr><td>MVP_score</td><td>The MVP score</td><td>Output</td></tr>

</table>
