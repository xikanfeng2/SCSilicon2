# Tutorials

## 1. Pre-requirements
* python3.6 or higher
* pandas>=0.23.4
* matplotlib>=3.0.2
* networkx>=3.2.1
* [wgsim](https://github.com/lh3/wgsim)


All python packages will be automatically installed when you install SCSilicon2 if these packages are not included in your python library.

To install wgsim, please refer to the README of [wgsim](https://github.com/lh3/wgsim). Please make sure the command 'wgsim' works in your command line.

## 2. Installation

### Creation of python virtual env
We recommend creating a virtual environment to run the scsilicon2(This step is optional!). You can use the following command to create a virtual python env:

```shell
# create a python virtual env in scsilicon2 folder
python -m venv scsilicon2

# activate the virtual env
source scsilicon2/bin/activate

# deactivate the virtual env
deactivate
```

### Installation with pip
To install with pip, run the following from a terminal:
```shell
pip install scsilicon2
```

### Installation from Github
To clone the repository and install manually, run the following from a terminal:
```shell
git clone https://github.com/xikanfeng2/SCSilicon2.git
cd SCSilicon2
python setup.py install
```

## 3. Quick start
The following code runs SCSilicon.

```
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

