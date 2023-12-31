# SCSilicon2

SCSilicon2: a single-cell genomics simulator which cost-effectively simulates single-cell genomics reads with haplotype-specific copy number annotation.

## 1. Pre-requirements
* python3.6 or higher
* pandas>=0.23.4
* matplotlib>=3.0.2
* networkx>=3.2.1
* [wgsim](https://github.com/lh3/wgsim)


All python packages will be automatically installed when you install SCSilicon2 if these packages are not included in your python library.

To install wgsim, please refer to the README of [wgsim](https://github.com/lh3/wgsim). Please make sure the command 'wgsim' works in your command line.

## 2. Installation

## Creation of python virtual env
We recommend creating a virtual environment to run the scsilicon2(This step is optional!). You can use the following command to create a virtual python env:

```Bash
# create a python virtual env in scsilicon2 folder
python -m venv scsilicon2

# activate the virtual env
source scsilicon2/bin/activate

# deactivate the virtual env
deactivate
```

### Installation with pip
To install with pip, run the following from a terminal:
```Bash
pip install scsilicon2
```

### Installation from Github
To clone the repository and install manually, run the following from a terminal:
```Bash
git clone https://github.com/xikanfeng2/SCSilicon2.git
cd SCSilicon2
python setup.py install
```

## 3. Quick start
The following code runs SCSilicon.

```Python
import scsilicon2 as scs

# create SCSilicon2 object: ref_genome and snp_file are required, and outdir, clone_no, and cell_no are optional.
simulator = scs.SCSilicon2(ref_genome='your reference fasta file here', snp_file='your snp list file here', outdir='your output directory here', clone_no=4, cell_no=10)

# simulate dataset
simulator.sim_dataset()
```

## 4. Input file required

1. **A reference genome file with fasta format.**  
Please refer to the example fasta file `example/input/chr22.fa`.
2. **A list of SNPs.**   
The SNPs in this list can be introduced in arbitrary positions of the genome. Please refer to the example snp list file `example/input/dbsnp.tsv`.

## 5. Output files of SCSilicon2
The output directory contains three subfolders: fastq folder, fasta folder and profile folder. The structure of one example output directory is listed as follows (the clone no is 3 and the cell no is 10 in this example):

```
output
 |-fastq
 | |-normal_r2.fq
 | |-clone2
 | | |-cell0_r1.fq
 | | |-cell0_r2.fq
 | |-normal_r1.fq
 | |-clone2_r2.fq
 | |-clone1_r1.fq
 | |-clone0_r2.fq
 | |-clone0
 | | |-cell2_r1.fq
 | | |-cell3_r2.fq
 | | |-cell2_r2.fq
 | | |-cell1_r2.fq
 | | |-cell1_r1.fq
 | | |-cell3_r1.fq
 | | |-cell0_r1.fq
 | | |-cell0_r2.fq
 | |-normal
 | | |-cell2_r1.fq
 | | |-cell2_r2.fq
 | | |-cell1_r2.fq
 | | |-cell1_r1.fq
 | | |-cell0_r1.fq
 | | |-cell0_r2.fq
 | |-clone2_r1.fq
 | |-clone1
 | | |-cell1_r2.fq
 | | |-cell1_r1.fq
 | | |-cell0_r1.fq
 | | |-cell0_r2.fq
 | |-clone0_r1.fq
 | |-clone1_r2.fq
 |-fasta
 | |-clone2.fasta
 | |-normal_paternal.fasta
 | |-clone2_paternal.fasta
 | |-clone0.fasta
 | |-clone1.fasta
 | |-clone0_paternal.fasta
 | |-normal.fasta
 | |-clone1_paternal.fasta
 | |-clone2_maternal.fasta
 | |-clone1_maternal.fasta
 | |-clone0_maternal.fasta
 | |-normal_maternal.fasta
 |-profile
 | |-changes.csv
 | |-tree.pdf
 | |-maternal_cnv_matrix.csv
 | |-paternal_cnv_matrix.csv
 | |-phases.csv
 | |-cnv_profile.csv
 | |-tree.newick
```

* `fasta folder`: stores all the fasta file for each clone.

* `fastq folder`: stores all the paired-reads with fastq format for each clone and each cell.

*  `profile folder`: stores all the profile file which is related to the simulation process. The detailed explanation of the format for each file in this folder is as follows.

    1. `changes.csv`: stores the evlution path for each clone. One example is listed below:

        |Parent|Child |Haplotype|Type|Segment                |Change|
        |------|------|---------|----|-----------------------|------|
        |normal|clone0|paternal |dup |chr22:500001-1000000   |1->3  |
        |normal|clone0|maternal |del |chr22:3500001-4000000  |1->0  |
        |normal|clone0|maternal |dup |chr22:4000001-4500000  |1->2  |
        |normal|clone0|maternal |dup |chr22:5000001-5500000  |1->2  |
        |normal|clone0|maternal |dup |chr22:8000001-8500000  |1->4  |
 

    2. `cnv_profile.csv`: stores the cnv ground truth for ech clone with maternal|paternal format. One example is listed below:

        |Chromosome|Start |End     |clone0|clone1                 |clone2|
        |----------|------|--------|------|-----------------------|------|
        |chr22     |1     |500000  |1&#124;1   |3&#124;1                    |3&#124;1   |
        |chr22     |500001|1000000 |1&#124;3   |1&#124;3                    |3&#124;5   |
        |chr22     |1000001|1500000 |1&#124;1   |3&#124;2                    |3&#124;2   |
        |chr22     |1500001|3000000 |1&#124;1   |1&#124;1                    |1&#124;1   |
        |chr22     |3000001|3500000 |1&#124;1   |3&#124;2                    |3&#124;2   |
 

    3. `maternal_cnv_matrix.csv` and `paternal_cnv_matrix.csv`: store the cnv matrix of each clone seperated by maternal haplotype and paternal haplotype. One example is listed below:

        |Index|clone0_maternal_cnvs|clone1_maternal_cnvs|clone2_maternal_cnvs|
        |------|--------------------|--------------------|--------------------|
        |chr22:1-500000|1                   |3                   |3                   |
        |chr22:500001-1000000|1                   |1                   |3                   |
        |chr22:1000001-1500000|1                   |3                   |3                   |
        |chr22:1500001-3000000|1                   |1                   |1                   |
        |chr22:3000001-3500000|1                   |3                   |3                   |
    
    4. `phases.csv`: stores the SNPs in maternal|paternal haplotype. One example is listed below:

        |chr22 |16578327|1&#124;0     |
        |------|--------|--------|
        |chr22 |17307398|1&#124;0     |
        |chr22 |18025718|1&#124;0     |
        |chr22 |21416314|0&#124;1     |
        |chr22 |22418251|1&#124;0     |

    5. `tree.newick` and `tree.pdf`: the cnv elution tree with newick format and pdf format.

    The example profile folder can be found in `data/profile` folder.

## 6. `SCSilicon2` object
All the general parameters for the SCSilicon2 simulation are stored in a `SCSilicon2` object. Let’s create a new one.

```Python
simulator = scs.SCSilicon2()
```

### 6.1 All parameters in `SCSilicon2` object

* `ref_genome`: str, required<br>
    The reference genome file path
        
* `snp_file`: str, required<br>
    The snp list file

* `outdir`: str, optional, default: './'<br>
    The output directory

* `clone_no`: int, optional, default: 1<br>
    The random clone number contained in evolution tree

* `cell_no`: int, optional, default: 2<br>
    The total cell number for this simultion dataset. Please make sure the `cell_no` is large than `clone_no`. At least one cell is geneated for nomal case.

* `max_cnv_tree_depth`: int, optional, default: 4<br>
    The maximum depth of random evolution tree

* `bin_len`: int, optional, default: 500000<br>
    The fixed bin length

* `HEHOratio`: float, optional, default: 0.5<br>
    Ratio of heterozygous SNPs

* `cnv_prob_cutoff`: float, optional, default: 0.8<br>
    The cutoff probability of a bin undergoing CNV, if random probability is larger than cutoff, CNV happens

* `clone_coverage`: float, optional, default: 30<br>
    The coverage for clone fastq file

* `cell_coverage`: float, optional, default: 0.5<br>
    The coverage for each cell in a clone

* `reads_len`: int, optional, default: 150<br>
    The reads length in fastq file

* `insertion_size`: int, optional, default: 350<br>
    The outer distance between the two ends

* `error_rate`: float, optional, default: 0.02<br>
    The base error rate

### 4.2 Getting and setting
If we want to look at the value of parameters, we can extract it using the `get_params` function:

```Python
simulator.get_params()

# console log: {'out_dir': './', 'ref': 'hg19', 'chrom': 'chr20', 'layout': 'PE', 'coverage': 5, 'isize': 260, 'threads': 10}
```

Alternatively, to give a parameter a new value we can use the `set_params` function:

```Python
simulator.set_params(clone_no=6, cell_no=20)
```

## 5. Simulating reads for SNPs using `SNPSimulator` object
Once we have a set of parameters we are happy with we can use `SNPSimulator` to simulate samples with SNPs in it. 
 ```Python
snp_simulator = scs.SNPSimulator()
snp_simulator.sim_samples(params)
```

### 5.1 All parameters in `SNPSimulator` object

* `cell_no`: int, optional, default: 1.<br>
    The cell number for this simulation

* snp_no : int, optional, default: 1000<br>
    The SNP number of each sample
        

For each sample, `SNPSimulator` will randomly select a total number of SNPs from dbSNP file and `snp_no` parameter can be used to control this total number.

### 5.2 Getting and setting
Similar to `SCSiliconParams`, `SNPSimulator` uses the functions `get_params` and `set_params` to get or set parameters.

### 5.3 Generating FASTAQ sample
`SNPSimulator` object uses the function `sim_samples` to generate FASTQ files for each sample. 
```Python
snp_simulator.sim_samples()
```
If you want to simulate `multiple` samples once, you can use the `cell_no` parameter to contorl this.
```Python
snp_simulator.set_params(cell_no=10)

# or set the parameter when creating the object
snp_simulator = scs.SNPSimulator(cell_no=10)

# generating reads
snp_simulator.sim_samples(params)
```
Above code will simulate 10 samples with FASTQ format once. 

### 5.4 Output files of `sim_samples` function
The `sim_samples` function will generate two output files for each sample in your output directory.

- `sample{1}-snps.txt`: the SNPs included in this sample. This file can be reagrded as the groud truth for SNP detection software. 
- `sample{1}.fq`: the reads data of this sample with FASTQ format.

`{1}` is the sample no., like sample1-snps.txt, sample2-snps.txt.


## 6. Simulating reads for CNVs using `CNVimulator` object
We can use `CNVimulator` to simulate samples with CNVs.

 ```Python
cnv_simulator = scs.CNVSimulator()
cnv_simulator.sim_samples(params)
```

### 6.1 All parameters in `CNVimulator` object

* `cell_no`: int, optional, default: 1.<br>
    The cell number for this simulation

* `bin_len`: int, optional, default: 500000.<br>
    The fixed bin length

* `seg_no`: int, optional, default: 10.<br>
    The segment number for each chromosome

* `cluster_no`: int, optional, default: 1.<br>
    The cell cluster number for multiple sample simulation

* `normal_frac`: float, optional, default: 0.4.<br>
    The fraction of normal cells

* `noise_frac`: float, optional, default: 0.1.<br>
    The noise fraction for cnv matrix

### 6.2 Getting and setting
Similar to `SCSiliconParams`, `CNVimulator` uses the functions `get_params` and `set_params` to get or set parameters.

### 6.3 Generating FASTAQ sample
`CNVimulator` object also uses the function `sim_samples` to generate FASTQ files for each sample. 
```Python
cnv_simulator.sim_samples(params)
```
The `seg_no` parameter can be used to control the segments in each chromosome.
```Python
cnv_simulator.set_params(seg_no=8)

# or set the parameter when creating the object
cnv_simulator = scs.SNPSimulator(seg_no=8)

# generating reads
cnv_simulator.sim_samples(params)
```
Above code will split each chromosome to 8 segments and this is useful for segmentation experiments of single cell CNV detection tools.

If you want to simulate `multiple` samples once, you can use the `cell_no` parameter to contorl this.
```Python
cnv_simulator.set_params(cell_no=10)

# or set the parameter when creating the object
cnv_simulator = scs.SNPSimulator(cell_no=10)

# generating reads
cnv_simulator.sim_samples(params)
```
Above code will simulate 10 samples with FASTQ format once.

For multiple-sample simulation, you can use the `cluster_no` parameter to seperate these samples to several clusters.
```Python
cnv_simulator.set_params(cluster_no=5)

# or set the parameter when creating the object
cnv_simulator = scs.SNPSimulator(cluster_no=10)

# generating reads
cnv_simulator.sim_samples(params)
```
### 6.4 Output files of `sim_samples` function
The `sim_samples` function will generate two output files for each sample in your output directory.

- `cnv.csv`: the CNV matrix with cells as rows and bins as columns. This file can be reagrded as the groud truth for CNV detection software. 
- `segments.csv`:  the segments information for each chromosome. This file can be reagrded as the groud truth for segmentation experiments.
- `clusters.csv`:  the clusters information for each sample. This file can be reagrded as the groud truth for cell cluster experiments.
- `sample{1}.fq`: the reads data of this sample with FASTQ format.

`{1}` is the sample no., like sample1.fq, sample2.fq.


### 6.5 Visualizing the CNV matrix
`CNVimulator` object has the funciton `visualize_cnv_matrix` to draw the heatmap graph for the cnv matrix.
```Python
cnv_simulator.visualize_cnv_matrix(out_prefix)
```
This function will save the heatmap with pdf format to the file named as `out_prefix.pdf`. One example of cnv heatmap graph is shown below:

![cnv heatmap](cnv-heatmap.png)


## 7. Simulating reads for SNVs using `SNVSimulator` object
Once we have a set of parameters we are happy with we can use `SNVSimulator` to simulate samples with SNVs in it. 
 ```Python
snv_simulator = scs.SNVSimulator()
snv_simulator.sim_samples(params)
```
### 7.1 All parameters in `SNVSimulator` object

* `cell_no`: int, optional, default: 1.<br>
    The cell number for this simulation

* `snv_no`: int, optional, default: 1000<br>
    The SNV number of each sample

        

## 8. Simulating reads for Indels using `IndelSimulator` object
Once we have a set of parameters we are happy with we can use `IndelSimulator` to simulate samples with Indels in it. 
 ```Python
indel_simulator = scs.IndelSimulator()
indel_simulator.sim_samples(params)
```
### 8.1 All parameters in `IndelSimulator` object

* `cell_no`: int, optional, default: 1.<br>
    The cell number for this simulation

* `in_no`: int, optional, default: 1000<br>
    The insertion number of each sample

* `del_no`: int, optional, default: 1000<br>
    The deletion number of each sample

## Cite us
todo

## Help
If you have any questions or require assistance using SCSilicon, please contact us with fxk@nwpu.edu.cn.