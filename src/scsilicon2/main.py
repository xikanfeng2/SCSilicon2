import os
import copy
import random
import pandas as pd
import numpy as np
import subprocess
import logging
from . import utils
from . import random_tree
from collections import deque
from glob import glob
from multiprocessing.pool import ThreadPool as Pool
from pathlib import Path

pd.options.mode.chained_assignment = None

class SCSilicon2:
    def __init__(self, ref_genome, snp_file=None, ignore_file=None, outdir='./', clone_no=1, cell_no=2, max_cnv_tree_depth=4, bin_len=500000, snp_ratio=0.000000333, wgsim_thread=1, HEHO_ratio=0.5, cnv_prob_cutoff=0.8, clone_coverage=15, cell_coverage=0.5, reads_len=150, insertion_size=350, error_rate=0.02, WGD_no=0, WCL_no=0, CNL_LOH_no=10, CNN_LOH_no=10, GOH_no=10, mirrored_cnv_no=10, barcodes_file=None, mode=0, bwa_thread=1, wgsim_path='wgsim', samtools_path='samtools', bwa_path='bwa', picard_path='picard.jar'):
        self.ref_genome = ref_genome
        self.snp_file = snp_file
        self.ignore_file = ignore_file
        self.outdir = outdir
        self.clone_no = clone_no
        self.cell_no = cell_no
        self.wgsim_thread = wgsim_thread
        self.bwa_thread = bwa_thread
        self.max_cnv_tree_depth = max_cnv_tree_depth
        self.bin_len = bin_len
        self.snp_ratio = snp_ratio
        self.HEHO_ratio = HEHO_ratio
        self.cnv_prob_cutoff = cnv_prob_cutoff
        self.clone_coverage = clone_coverage
        self.cell_coverage = cell_coverage
        self.reads_len = reads_len
        self.insertion_size = insertion_size
        self.error_rate = error_rate
        self.WGD_no = WGD_no
        self.WCL_no = WCL_no
        self.CNL_LOH_no = CNL_LOH_no
        self.CNN_LOH_no = CNN_LOH_no
        self.GOH_no = GOH_no
        self.mirrored_cnv_no = mirrored_cnv_no
        self.mode = mode
        if barcodes_file == None:
            barcodes_file = os.path.join(utils.root_path(), 'data/barcodes.txt')
        else:
            barcodes_file = barcodes_file
        self.wgsim_path = wgsim_path
        self.samtools_path = samtools_path
        self.bwa_path = bwa_path
        self.picard_path = picard_path
        self.chrom_sizes = {}
        self.ignore_list = []
        self._check_params()
        self.samples = dict.fromkeys(['cell' + str(i+1) for i in range(self.cell_no)])
        for sample in self.samples:
            self.samples[sample] = {}
        
        # config log
        log_file = os.path.join(self.outdir, 'log.txt')
        logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.DEBUG, handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ])


    def _check_params(self):
        """Check SCSilicon parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as threads='10.5', would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        utils.check_exist(ref_genome=self.ref_genome)
        # utils.check_exist(snp_file=self.snp_file)
        utils.check_int(clone_no=self.clone_no)
        utils.check_positive(clone_no=self.clone_no)
        utils.check_int(cell_no=self.cell_no)
        utils.check_positive(cell_no=self.cell_no)
        if self.cell_no <= self.clone_no:
            raise ValueError(
                "cell_no should larger than clone_no.")
        utils.check_int(max_cnv_tree_depth=self.max_cnv_tree_depth)
        utils.check_positive(max_cnv_tree_depth=self.max_cnv_tree_depth)
        utils.check_int(bin_len=self.bin_len)
        utils.check_positive(bin_len=self.bin_len)
        utils.check_between(0,1,HEHO_ratio=self.snp_ratio)
        utils.check_between(0,1,HEHO_ratio=self.HEHO_ratio)
        utils.check_between(0,1,cnv_prob_cutoff=self.cnv_prob_cutoff)
        utils.check_positive(clone_coverage=self.clone_coverage)
        utils.check_positive(cell_coverage=self.cell_coverage)
        utils.check_int(reads_len=self.reads_len)
        utils.check_positive(reads_len=self.reads_len)
        utils.check_int(reads_len=self.wgsim_thread)
        utils.check_positive(reads_len=self.wgsim_thread)
        utils.check_int(insertion_size=self.insertion_size)
        utils.check_positive(insertion_size=self.insertion_size)
        utils.check_between(0,1,error_rate=self.error_rate)
        utils.check_int(WGD_no=self.WGD_no)
        utils.check_lt_zero(WGD_no=self.WGD_no)
        utils.check_int(WCL_no=self.WCL_no)
        utils.check_lt_zero(WCL_no=self.WCL_no)
        utils.check_int(CNL_LOH_no=self.CNL_LOH_no)
        utils.check_lt_zero(CNL_LOH_no=self.CNL_LOH_no)
        utils.check_int(CNN_LOH_no=self.CNN_LOH_no)
        utils.check_lt_zero(CNN_LOH_no=self.CNN_LOH_no)
        utils.check_int(GOH_no=self.GOH_no)
        utils.check_lt_zero(GOH_no=self.GOH_no)
        utils.check_int(mirrored_cnv_no=self.mirrored_cnv_no)
        utils.check_lt_zero(mirrored_cnv_no=self.mirrored_cnv_no)
    def set_params(self, **params):
        """Set the parameters of SCSilicon2.

        Any parameters not given as named arguments will be left at their
        current value.

        Parameters
        ----------

        ref_genome: str, required
            The reference genome file path
        
        snp_file: str, required
            The snp list file

        outdir: str, optional, default: './'
            The output directory
        
        clone_no: int, optional, default: 1
            The random clone number contained in evolution tree
        
        cell_no: int, optional, default: 2
            The total cell number for this simultion dataset
        
        max_cnv_tree_depth: int, optional, default: 4
            The maximum depth of random evolution tree
        
        bin_len: int, optional, default: 500000
            The fixed bin length
        
        HEHO_ratio: float, optional, default: 0.5
            Ratio of heterozygous SNPs

        WGD_no: int, optional, default: 0

        
        cnv_prob_cutoff: float, optional, default: 0.8
            The cutoff probability of a bin undergoing CNV, if random probability is larger than cutoff, CNV happens
        
        clone_coverage: float, optional, default: 30
            The coverage for clone fastq file

        cell_coverage: float, optional, default: 0.5
            The coverage for each cell in a clone
        
        reads_len: int, optional, default: 150
            The reads length in fastq file
        
        insertion_size: int, optional, default: 350
            The outer distance between the two ends
        
        error_rate: float, optional, default: 0.02
            The base error rate

        Returns
        -------
        self
        """

        # parameters
        if 'ref_genome' in params and params['ref_genome'] != self.ref_genome:
            self.ref_genome = params['ref_genome']
            del params['ref_genome']
        if 'snp_file' in params and params['snp_file'] != self.snp_file:
            self.snp_file = params['snp_file']
            del params['snp_file']
        if 'ignore_file' in params and params['ignore_file'] != self.ignore_file:
            self.ignore_file = params['ignore_file']
            del params['ignore_file']
        if 'outdir' in params and params['outdir'] != self.outdir:
            self.outdir = params['outdir']
            del params['outdir']
        if 'clone_no' in params and params['clone_no'] != self.clone_no:
            self.clone_no = params['clone_no']
            del params['clone_no']
        if 'cell_no' in params and params['cell_no'] != self.cell_no:
            self.cell_no = params['cell_no']
            del params['cell_no']
        if 'max_cnv_tree_depth' in params and params['max_cnv_tree_depth'] != self.max_cnv_tree_depth:
            self.max_cnv_tree_depth = params['max_cnv_tree_depth']
            del params['max_cnv_tree_depth']
        if 'bin_len' in params and params['bin_len'] != self.bin_len:
            self.bin_len = params['bin_len']
            del params['bin_len']
        if 'wgsim_thread' in params and params['wgsim_thread'] != self.wgsim_thread:
            self.wgsim_thread = params['wgsim_thread']
            del params['wgsim_thread']
        if 'HEHO_ratio' in params and params['HEHO_ratio'] != self.HEHO_ratio:
            self.HEHO_ratio = params['HEHO_ratio']
            del params['HEHO_ratio']
        if 'cnv_prob_cutoff' in params and params['cnv_prob_cutoff'] != self.cnv_prob_cutoff:
            self.cnv_prob_cutoff = params['cnv_prob_cutoff']
            del params['cnv_prob_cutoff']
        if 'clone_coverage' in params and params['clone_coverage'] != self.clone_coverage:
            self.clone_coverage = params['clone_coverage']
            del params['clone_coverage']
        if 'cell_coverage' in params and params['cell_coverage'] != self.cell_coverage:
            self.cell_coverage = params['cell_coverage']
            del params['cell_coverage']
        if 'reads_len' in params and params['reads_len'] != self.reads_len:
            self.reads_len = params['reads_len']
            del params['reads_len']
        if 'insertion_size' in params and params['insertion_size'] != self.insertion_size:
            self.insertion_size = params['insertion_size']
            del params['insertion_size']
        if 'error_rate' in params and params['error_rate'] != self.error_rate:
            self.error_rate = params['error_rate']
            del params['error_rate']
        if 'WGD_no' in params and params['WGD_no'] != self.WGD_no:
            self.WGD_no = params['WGD_no']
            del params['WGD_no']
        if 'WCL_no' in params and params['WCL_no'] != self.WCL_no:
            self.WCL_no = params['WCL_no']
            del params['WCL_no']
        if 'CNL_LOH_no' in params and params['CNL_LOH_no'] != self.CNL_LOH_no:
            self.CNL_LOH_no = params['CNL_LOH_no']
            del params['CNL_LOH_no']
        if 'CNN_LOH_no' in params and params['CNN_LOH_no'] != self.CNN_LOH_no:
            self.CNN_LOH_no = params['CNN_LOH_no']
            del params['CNN_LOH_no']
        if 'GOH_no' in params and params['GOH_no'] != self.GOH_no:
            self.GOH_no = params['GOH_no']
            del params['GOH_no']
        if 'mirrored_cnv_no' in params and params['mirrored_cnv_no'] != self.mirrored_cnv_no:
            self.mirrored_cnv_no = params['mirrored_cnv_no']
            del params['mirrored_cnv_no']
        self._check_params()
        self.get_params()
        return self

    def get_params(self):
        print(vars(self))

    def _get_chrom_sizes(self):
        if self.ignore_file:
            self.ignore_list = utils.parseIgnoreList(self.ignore_file)

        with open(self.ref_genome, 'r') as refinput:
            chrom = None
            for line in refinput:
                line = line.strip()
                if line.startswith('>'):
                    chrom = line.strip()[1:].split()[0]
                    if chrom not in self.ignore_list:
                        self.chrom_sizes[chrom] = 0
                else:
                    if chrom in self.ignore_list:
                        continue
                    linelen = len(line.strip())
                    self.chrom_sizes[chrom] += linelen

    def _buildGenome(self, maternalFasta, paternalFasta, phaselist):
        if self.snp_file == None:
            allsnps = utils.randomSNPList(self.chrom_sizes, self.snp_ratio)
        else:
            allsnps = utils.parseSNPList(self.snp_file)
        print(allsnps)
        phases = {}
        # m_genome = {}
        # p_genome = {}
        with open(self.ref_genome, 'r') as refinput:
            with open(maternalFasta, 'w') as out1:
                with open(paternalFasta, 'w') as out2:
                    chrom = None
                    snps = None
                    
                    for line in refinput:
                        line = line.strip()
                        if line.startswith('>'):
                            if chrom and chrom not in self.ignore_list:
                                out1.write('\n')
                                out2.write('\n')
                            chrom = line.strip()[1:].split()[0]
                            if chrom in self.ignore_list:
                                continue
                            out1.write(line+'\n')
                            out2.write(line+'\n')
                            # m_genome[chrom] = ''
                            # p_genome[chrom] = ''
                            snps = allsnps[chrom]
                            snppos = sorted(snps.keys())
                            currentpos = 0 
                            currentsnppos = snppos.pop(0)
                            allele1 = snps[currentsnppos][0]
                            allele2 = snps[currentsnppos][1]
                        else:
                            if chrom in self.ignore_list:
                                continue
                            linelen = len(line.strip())
                            if int(currentsnppos) > currentpos and int(currentsnppos) <= currentpos + linelen:
                                while int(currentsnppos) > currentpos and int(currentsnppos) <= currentpos + linelen:
                                    sindex = int(currentsnppos)-currentpos-1
                                    a = line[sindex]
                                    mline = line
                                    pline = line
                                    if random.random() < self.HEHO_ratio: #Heterozygous
                                        if random.random() < 0.5:
                                            a1 = allele1.lower() if a.islower() else allele1.upper()
                                            a2 = allele2.lower() if a.islower() else allele2.upper()
                                            phases[(chrom, currentsnppos)] = a1.upper() + ',' + a2.upper() + ',0|1'
                                            mline = mline[:sindex]+a1+mline[sindex+1:]
                                            pline = pline[:sindex]+a2+pline[sindex+1:]
                                        else:
                                            a1 = allele2.lower() if a.islower() else allele2.upper()
                                            a2 = allele1.lower() if a.islower() else allele1.upper()
                                            phases[(chrom, currentsnppos)] = a2.upper() + ',' + a1.upper() + ',1|0'
                                            mline = mline[:sindex]+a1+mline[sindex+1:]
                                            pline = pline[:sindex]+a2+pline[sindex+1:]
                                    else: #Homozygous
                                        a1 = allele1.lower() if a.islower() else allele1.upper()
                                        mline = mline[:sindex]+a1+mline[sindex+1:]
                                        pline = pline[:sindex]+a1+mline[sindex+1:]
                                    if snppos:
                                        currentsnppos = snppos.pop(0)
                                        allele1 = snps[currentsnppos][0]
                                        allele2 = snps[currentsnppos][1]
                                    else:
                                        break
                                # m_genome[chrom] += mline.strip()
                                # p_genome[chrom] += pline.strip()
                                out1.write(mline)
                                out2.write(pline)
                            else:
                                # m_genome[chrom] += line.strip()
                                # p_genome[chrom] += line.strip()
                                out1.write(line)
                                out2.write(line)
                        currentpos += len(line)
                    out1.write('\n')
                    out2.write('\n')
        with open(phaselist, 'w') as output:
            for g in sorted(phases.keys(), key=(lambda x : (int(''.join([l for l in x[0] if l.isdigit()])), x[1]))):
                output.write('{},{},{}\n'.format(g[0], g[1], phases[g]))

    def _split_chr_to_bins(self, chrom):
        """Split chromosomes to fixed-lenght bins

        Parameters
        ----------
        bin_len : int
            fixed-bin-length

        Returns
        -------
        ref: Dataframe of pandas
        """
        ref = pd.DataFrame(
            columns=['Chromosome', 'Start', 'End'])
        bin_len = self.bin_len
        if chrom != 'all':
            chrom_size = self.chrom_sizes[chrom]
            start = 1
            end = bin_len
            count = 1
            while(start < chrom_size):
                ref = pd.concat([ref, pd.DataFrame([{
                    'Chromosome': chrom,
                    'Start': start,
                    'End': min(end, chrom_size),
                }])], ignore_index=True)
                count += 1
                start = end + 1
                end = bin_len * count     
        else:
            for chrom, chrom_size in self.chrom_sizes.items():
                start = 1
                end = bin_len
                count = 1
                while(start < chrom_size):
                    ref = pd.concat([ref, pd.DataFrame([{
                        'Chromosome': chrom,
                        'Start': start,
                        'End': min(end, chrom_size),
                    }])], ignore_index=True)
                    count += 1
                    start = end + 1
                    end = bin_len * count            
        return ref

    def _generate_cnv_profile_for_each_clone(self, root, ref, m_fasta, p_fasta):
        cutoff = self.cnv_prob_cutoff
        changes = []
        all_chroms = list(np.unique(ref['Chromosome']))
        if self.WGD_no + self.WCL_no > len(all_chroms):
            raise Exception("The sum of WGD_no and WCL_no should be less or equal to the total number of chromosomes!")

        # store maternal and paternal genome to dict
        maternal_genome = {}
        paternal_genome = {}
        with open(m_fasta, 'r') as input:
            chrom = None
            for line in input:
                line = line.strip()
                if line.startswith('>'):
                    chrom = line.strip()[1:].split()[0]
                    maternal_genome[chrom] = ''
                else:
                    maternal_genome[chrom] += line
        with open(p_fasta, 'r') as input:
            chrom = None
            for line in input:
                line = line.strip()
                if line.startswith('>'):
                    chrom = line.strip()[1:].split()[0]
                    paternal_genome[chrom] = ''
                else:
                    paternal_genome[chrom] += line
        
        # add nromal clone to cnv matrix
        for i in range(ref.shape[0]):
            root.maternal_cnvs.append(1)
            root.paternal_cnvs.append(1)
        ref[root.name+'_maternal_cnvs'] = root.maternal_cnvs
        ref[root.name+'_paternal_cnvs'] = root.paternal_cnvs
        
        # add the children of normal clone to queue
        queue = deque(root.children)

        while  queue:
            clone = queue.popleft()
            if clone.depth == 1: # children of normal clone
                mirrored_cnv_flag = False

                wgd_chroms = []
                wcl_chroms = []
                
                # select WGD and WCL chromosomes
                random_chroms = random.sample(all_chroms, self.WGD_no+self.WCL_no)
                wgd_chroms = random_chroms[:self.WGD_no]
                wcl_chroms = random_chroms[self.WCL_no:]
                wgd_cnvs = dict.fromkeys(wgd_chroms) # store the cnv number for each wgd chrom
                wcl_cnvs = dict.fromkeys(wcl_chroms) # store the cnv number for each wgd chrom

                # select the position for CNL_LOH, CNN_LOH, GOH and mirrored cnv
                random_bins = random.sample(range(0, ref.shape[0]), self.CNL_LOH_no + self.CNN_LOH_no + self.GOH_no + self.mirrored_cnv_no)
                cnl_loh_bins = random_bins[:self.CNL_LOH_no]
                cnn_loh_bins = random_bins[self.CNL_LOH_no:self.CNL_LOH_no + self.CNN_LOH_no]
                goh_bins = random_bins[self.CNL_LOH_no + self.CNN_LOH_no:self.CNL_LOH_no + self.CNN_LOH_no + self.GOH_no]
                mirrored_cnv_bins = random_bins[self.CNL_LOH_no + self.CNN_LOH_no + self.GOH_no:]
                
                for i in range(ref.shape[0]):
                    # if flag is ture, the previous bin has been process as mirrored cnv bin and skip it.
                    if mirrored_cnv_flag:
                        mirrored_cnv_flag = False
                        continue
                    
                    current_chrom = ref['Chromosome'][i]
                    start = ref['Start'][i]
                    end = ref['End'][i]
                    m_sequence = maternal_genome[current_chrom][start-1:end]
                    p_sequence = paternal_genome[current_chrom][start-1:end]

                    # handle WGD
                    if current_chrom in wgd_chroms:
                        if not wgd_cnvs[current_chrom]:
                            m_cnv = utils.random_WGD()
                            p_cnv = utils.random_WGD()
                            wgd_cnvs[current_chrom] = (m_cnv, p_cnv)
                        else:
                            m_cnv = wgd_cnvs[current_chrom][0]
                            p_cnv = wgd_cnvs[current_chrom][1]

                        # 1. WGD in maternal and paternal 2. WGD in maternal 3. WGD in paternal
                        random_prob = random.random()
                        if random_prob < 1/3:
                            changes.append(['normal',clone.name,'maternal','WGD',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(m_cnv)])
                            changes.append(['normal',clone.name,'paternal','WGD',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(p_cnv)])
                            clone.maternal_cnvs.append(m_cnv)
                            clone.paternal_cnvs.append(p_cnv)
                        elif random_prob > 2/3:
                            changes.append(['normal',clone.name,'maternal','WGD',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(m_cnv)])
                            clone.maternal_cnvs.append(m_cnv)
                            clone.paternal_cnvs.append(1)
                        else:
                            changes.append(['normal',clone.name,'paternal','WGD',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(p_cnv)])
                            clone.maternal_cnvs.append(1)
                            clone.paternal_cnvs.append(p_cnv)
                        clone.changes.append('WGD')
                        continue
                    
                    # handle WCL
                    if current_chrom in wcl_chroms:
                        if not wcl_cnvs[current_chrom]:
                            # 1. WCL in maternal and paternal 2. WCL in maternal 3. WCL in paternal
                            random_prob = random.random()
                            if random_prob < 1/3:
                                wcl_cnvs[current_chrom] = (True, True)
                            elif random_prob > 2/3:
                                wcl_cnvs[current_chrom] = (True, False)
                            else:
                                wcl_cnvs[current_chrom] = (False, True)

                        if wcl_cnvs[current_chrom][0]:
                            changes.append(['normal',clone.name,'maternal','WCL',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(0)])
                            clone.maternal_cnvs.append(0)
                        else:
                            clone.maternal_cnvs.append(1)
                        if wcl_cnvs[current_chrom][1]:
                            changes.append(['normal',clone.name,'paternal','WCL',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(0)])
                            clone.paternal_cnvs.append(0)
                        else:
                            clone.paternal_cnvs.append(1)
                        clone.changes.append('WCL')
                        continue
                    
                    # handle CNL_LOH: 1:0 or 0:1
                    if i in cnl_loh_bins:
                        # check heterozygosity
                        if m_sequence != p_sequence:
                            cnl_cnv = utils.random_CNL()

                            if random.random() < 0.5: # m:p = 1:0
                                if cnl_cnv != 1:
                                    changes.append(['normal',clone.name,'maternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnl_cnv)])
                                changes.append(['normal',clone.name,'paternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(0)])
                                clone.maternal_cnvs.append(cnl_cnv)
                                clone.paternal_cnvs.append(0)
                            else: # m:p = 0:1
                                if cnl_cnv != 1:
                                    changes.append(['normal',clone.name,'paternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnl_cnv)])
                                changes.append(['normal',clone.name,'maternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(0)])
                                clone.maternal_cnvs.append(0)
                                clone.paternal_cnvs.append(cnl_cnv)
                            clone.changes.append('CNL_LOH')
                            continue
                        else:
                            clone.maternal_cnvs.append(1)
                            clone.paternal_cnvs.append(1)
                            clone.changes.append('NONE')
                            continue
                    
                    # handle CNN_LOH: 2:0 or 0:2
                    if i in cnn_loh_bins:
                        # check heterozygosity
                        if m_sequence != p_sequence:
                            cnn_cnv = utils.random_mirrored_cnv()

                            if random.random() < 0.5: # m:p = 2:0
                                changes.append(['normal',clone.name,'maternal','CNN_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnn_cnv)])
                                changes.append(['normal',clone.name,'paternal','CNN_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(0)])
                                clone.maternal_cnvs.append(cnn_cnv)
                                clone.paternal_cnvs.append(0)
                            else: # m:p = 0:1
                                changes.append(['normal',clone.name,'maternal','CNN_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(0)])
                                changes.append(['normal',clone.name,'paternal','CNN_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnn_cnv)])
                                clone.maternal_cnvs.append(0)
                                clone.paternal_cnvs.append(cnn_cnv)
                            clone.changes.append('CNN_LOH')
                            continue
                        else:
                            clone.maternal_cnvs.append(1)
                            clone.paternal_cnvs.append(1)
                            clone.changes.append('NONE')
                            continue
                    
                    # handle GOH
                    if i in goh_bins:
                        # check heterozygosity
                        if m_sequence == p_sequence:
                            m_cnv = utils.random_WGD()
                            p_cnv = utils.random_WGD()
                            changes.append(['normal',clone.name,'maternal','GOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(m_cnv)])
                            changes.append(['normal',clone.name,'paternal','GOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(p_cnv)])
                            clone.maternal_cnvs.append(m_cnv)
                            clone.paternal_cnvs.append(p_cnv)
                            clone.changes.append('GOH')
                            continue
                        else:
                            clone.maternal_cnvs.append(1)
                            clone.paternal_cnvs.append(1)
                            clone.changes.append('NONE')
                            continue
                    
                    # handle mirrored cnv
                    if i in mirrored_cnv_bins:
                        # make sure the next bin located in same chromosome
                        if i+1 >= ref.shape[0] or ref['Chromosome'][i] != ref['Chromosome'][i+1]:
                            clone.maternal_cnvs.append(1)
                            clone.paternal_cnvs.append(1)
                            clone.changes.append('NONE')
                            continue

                        # generate mirrored cnv number
                        total_cnv = utils.random_mirrored_cnv()
                        cnv1 = random.randint(0,total_cnv)
                        while cnv1 == total_cnv/2:
                            cnv1 = random.randint(0, total_cnv)
                        cnv2 = total_cnv - cnv1
                        if random.random() < 0.5: # m:p = cnv1:cnv2
                            changes.append(['normal',clone.name,'maternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnv1)])
                            changes.append(['normal',clone.name,'maternal','mirrored cnv',ref['Chromosome'][i+1]+':'+str(ref['Start'][i+1])+'-'+str(ref['End'][i+1]),'1->'+str(cnv2)])
                            changes.append(['normal',clone.name,'paternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnv2)])
                            changes.append(['normal',clone.name,'paternal','mirrored cnv',ref['Chromosome'][i+1]+':'+str(ref['Start'][i+1])+'-'+str(ref['End'][i+1]),'1->'+str(cnv1)])
                            clone.maternal_cnvs.append(cnv1)
                            clone.paternal_cnvs.append(cnv2)
                            clone.maternal_cnvs.append(cnv2)
                            clone.paternal_cnvs.append(cnv1)
                        else: # m:p = cnv2:cnv1
                            changes.append(['normal',clone.name,'maternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnv2)])
                            changes.append(['normal',clone.name,'maternal','mirrored cnv',ref['Chromosome'][i+1]+':'+str(ref['Start'][i+1])+'-'+str(ref['End'][i+1]),'1->'+str(cnv1)])
                            changes.append(['normal',clone.name,'paternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnv1)])
                            changes.append(['normal',clone.name,'paternal','mirrored cnv',ref['Chromosome'][i+1]+':'+str(ref['Start'][i+1])+'-'+str(ref['End'][i+1]),'1->'+str(cnv2)])
                            clone.maternal_cnvs.append(cnv2)
                            clone.paternal_cnvs.append(cnv1)
                            clone.maternal_cnvs.append(cnv1)
                            clone.paternal_cnvs.append(cnv2)
                        mirrored_cnv_flag = True
                        clone.changes.append('mirrored cnv')
                        clone.changes.append('mirrored cnv')
                        continue
                    
                    # generate random cnv
                    if random.random() > cutoff: # 20% cnv
                        m_cnv = utils.random_cnv()
                        p_cnv = utils.random_cnv()
                        

                        # check whether is CNL_LOH
                        if (m_sequence != p_sequence) and ((m_cnv == 0 and p_cnv !=0) or (m_cnv != 0 and p_cnv ==0)):
                            changes.append(['normal',clone.name,'maternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(m_cnv)])
                            changes.append(['normal',clone.name,'paternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(p_cnv)])
                            clone.maternal_cnvs.append(m_cnv)
                            clone.paternal_cnvs.append(p_cnv)
                            clone.changes.append('CNL_LOH')
                            continue
                        
                        # check mirrored CNV
                        if clone.changes and clone.changes[-1] in ['REGULAR', 'NONE']:
                            if m_cnv != p_cnv and m_cnv == clone.paternal_cnvs[-1] and clone.maternal_cnvs[-1] == p_cnv:
                                if clone.changes[-1] == 'REGULAR':
                                    clone.changes.pop()
                                    changes.pop()
                                    changes.pop()
                                changes.append(['normal',clone.name,'maternal','mirrored cnv',ref['Chromosome'][i-1]+':'+str(ref['Start'][i-1])+'-'+str(ref['End'][i-1]),'1->'+str(clone.maternal_cnvs[-1])])
                                changes.append(['normal',clone.name,'maternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(m_cnv)])
                                changes.append(['normal',clone.name,'paternal','mirrored cnv',ref['Chromosome'][i-1]+':'+str(ref['Start'][i-1])+'-'+str(ref['End'][i-1]),'1->'+str(clone.paternal_cnvs[-1])])
                                changes.append(['normal',clone.name,'paternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(p_cnv)])
                                clone.maternal_cnvs.append(m_cnv)
                                clone.paternal_cnvs.append(p_cnv)
                                clone.changes.append('mirrored cnv')
                                clone.changes.append('mirrored cnv')
                                continue
                        
                        # normal case
                        if m_cnv == 0:
                            mtype = 'del'
                        else:
                            mtype = 'dup'
                        changes.append(['normal',clone.name,'maternal',mtype,ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(m_cnv)])
                        if p_cnv == 0:
                            ptype = 'del'
                        else:
                            ptype = 'dup'
                        changes.append(['normal',clone.name,'paternal',ptype,ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(p_cnv)])

                        clone.maternal_cnvs.append(m_cnv)
                        clone.paternal_cnvs.append(p_cnv)
                        clone.changes.append('REGULAR')
                    else: # normal
                        clone.maternal_cnvs.append(1)
                        clone.paternal_cnvs.append(1)
                        clone.changes.append('NONE')
            else:
                # first inherit from parent
                clone.maternal_cnvs = copy.deepcopy(clone.parent.maternal_cnvs)
                clone.paternal_cnvs = copy.deepcopy(clone.parent.paternal_cnvs)
                clone.changes = copy.deepcopy(clone.parent.changes)

                 # select the position for CNL_LOH, CNN_LOH, GOH and mirrored cnv
                random_bins = random.sample(range(0, ref.shape[0]), self.CNL_LOH_no + self.CNN_LOH_no + self.GOH_no + self.mirrored_cnv_no)
                cnl_loh_bins = random_bins[:self.CNL_LOH_no]
                cnn_loh_bins = random_bins[self.CNL_LOH_no:self.CNL_LOH_no + self.CNN_LOH_no]
                goh_bins = random_bins[self.CNL_LOH_no + self.CNN_LOH_no:self.CNL_LOH_no + self.CNN_LOH_no + self.GOH_no]
                mirrored_cnv_bins = random_bins[self.CNL_LOH_no + self.CNN_LOH_no + self.GOH_no:]
                mirrored_cnv_flag = False

                for i, item in enumerate(clone.paternal_cnvs):
                    # if flag is ture, the previous bin has been process as mirrored cnv bin and skip it.
                    if mirrored_cnv_flag:
                        mirrored_cnv_flag = False
                        continue
                    
                    if clone.parent.changes[i] not in ['REGULAR', 'NONE']:
                        continue

                    current_chrom = ref['Chromosome'][i]
                    start = ref['Start'][i]
                    end = ref['End'][i]
                    m_sequence = maternal_genome[current_chrom][start-1:end]
                    p_sequence = paternal_genome[current_chrom][start-1:end]

                    m_parent_cnv = clone.maternal_cnvs[i]
                    p_parent_cnv = clone.maternal_cnvs[i]

                    if clone.parent.changes[i] == 'NONE':
                        # handle CNL_LOH: 1:0 or 0:1
                        if i in cnl_loh_bins:
                            # check heterozygosity
                            if m_sequence != p_sequence:
                                cnl_cnv = utils.random_CNL()

                                if random.random() < 0.5: # m:p = 1:0
                                    if cnl_cnv != 1:
                                        changes.append([clone.parent.name,clone.name,'maternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnl_cnv)])
                                    changes.append([clone.parent.name,clone.name,'paternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(0)])
                                    clone.maternal_cnvs[i] = cnl_cnv
                                    clone.paternal_cnvs[i] = 0
                                else: # m:p = 0:1
                                    if cnl_cnv != 1:
                                        changes.append([clone.parent.name,clone.name,'paternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnl_cnv)])
                                    changes.append([clone.parent.name,clone.name,'maternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(0)])
                                    clone.maternal_cnvs[i] = 0
                                    clone.paternal_cnvs[i] = cnl_cnv
                                clone.changes[i] = 'CNL_LOH'
                                continue
                            else:
                                continue
                        
                        # handle CNN_LOH: 2:0 or 0:2
                        if i in cnn_loh_bins:
                            # check heterozygosity
                            if m_sequence != p_sequence:
                                cnn_cnv = utils.random_mirrored_cnv()

                                if random.random() < 0.5: # m:p = 2:0
                                    changes.append([clone.parent.name,clone.name,'maternal','CNN_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnn_cnv)])
                                    changes.append([clone.parent.name,clone.name,'paternal','CNN_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(0)])
                                    clone.maternal_cnvs[i] = cnn_cnv
                                    clone.paternal_cnvs[i] = 0
                                else: # m:p = 0:1
                                    changes.append([clone.parent.name,clone.name,'maternal','CNN_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(0)])
                                    changes.append([clone.parent.name,clone.name,'paternal','CNN_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnn_cnv)])
                                    clone.maternal_cnvs[i] = 0
                                    clone.paternal_cnvs[i] = cnn_cnv
                                clone.changes[i] = 'CNN_LOH'
                                continue
                            else:
                                continue
                    
                        # handle GOH
                        if i in goh_bins:
                            # check heterozygosity
                            if m_sequence == p_sequence:
                                m_cnv = utils.random_WGD()
                                p_cnv = utils.random_WGD()
                                changes.append([clone.parent.name,clone.name,'maternal','GOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(m_cnv)])
                                changes.append([clone.parent.name,clone.name,'paternal','GOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(p_cnv)])
                                clone.maternal_cnvs[i] = m_cnv
                                clone.paternal_cnvs[i] = p_cnv
                                clone.changes[i] = 'GOH'
                                continue
                            else:
                                continue
                    
                        # handle mirrored cnv
                        if i in mirrored_cnv_bins:
                            # make sure the next bin located in same chromosome
                            if i+1 >= len(clone.paternal_cnvs) or ref['Chromosome'][i] != ref['Chromosome'][i+1] or clone.parent.changes[i+1] != 'NONE':
                                continue

                            # generate mirrored cnv number
                            total_cnv = utils.random_mirrored_cnv()
                            cnv1 = random.randint(0,total_cnv)
                            while cnv1 == total_cnv/2:
                                cnv1 = random.randint(0, total_cnv)
                            cnv2 = total_cnv - cnv1
                            if random.random() < 0.5: # m:p = cnv1:cnv2
                                changes.append([clone.parent.name,clone.name,'maternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnv1)])
                                changes.append([clone.parent.name,clone.name,'maternal','mirrored cnv',ref['Chromosome'][i+1]+':'+str(ref['Start'][i+1])+'-'+str(ref['End'][i+1]),'1->'+str(cnv2)])
                                changes.append([clone.parent.name,clone.name,'paternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnv2)])
                                changes.append([clone.parent.name,clone.name,'paternal','mirrored cnv',ref['Chromosome'][i+1]+':'+str(ref['Start'][i+1])+'-'+str(ref['End'][i+1]),'1->'+str(cnv1)])
                                clone.maternal_cnvs[i] = cnv1
                                clone.paternal_cnvs[i] = cnv2
                                clone.maternal_cnvs[i+1] = cnv2
                                clone.paternal_cnvs[i+1] = cnv1
                            else: # m:p = cnv2:cnv1
                                changes.append([clone.parent.name,clone.name,'maternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnv2)])
                                changes.append([clone.parent.name,clone.name,'maternal','mirrored cnv',ref['Chromosome'][i+1]+':'+str(ref['Start'][i+1])+'-'+str(ref['End'][i+1]),'1->'+str(cnv1)])
                                changes.append([clone.parent.name,clone.name,'paternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(cnv1)])
                                changes.append([clone.parent.name,clone.name,'paternal','mirrored cnv',ref['Chromosome'][i+1]+':'+str(ref['Start'][i+1])+'-'+str(ref['End'][i+1]),'1->'+str(cnv2)])
                                clone.maternal_cnvs[i] = cnv2
                                clone.paternal_cnvs[i] = cnv1
                                clone.maternal_cnvs[i+1] = cnv1
                                clone.paternal_cnvs[i+1] = cnv2
                            mirrored_cnv_flag = True
                            clone.changes[i] = 'mirrored cnv'
                            clone.changes[i+1] = 'mirrored cnv'
                            continue
                    
                    # regular situation
                    if random.random() > cutoff: # 20% cnv
                        m_cnv = utils.random_cnv()
                        p_cnv = utils.random_cnv()
                        if m_parent_cnv == 0:
                            m_cnv = 0
                        else:
                            m_cnv = random.randint(m_parent_cnv, max(5, m_parent_cnv+1))
                        
                        if p_parent_cnv == 0:
                            p_cnv = 0
                        else:
                            p_cnv = random.randint(p_parent_cnv, max(5, p_parent_cnv+1))
                        
                        # check whether is CNL_LOH
                        if (m_sequence != p_sequence) and ((m_cnv == 0 and p_cnv !=0) or (m_cnv != 0 and p_cnv ==0)):
                            changes.append([clone.parent.name,clone.name,'maternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),str(m_parent_cnv)+'->'+str(m_cnv)])
                            changes.append([clone.parent.name,clone.name,'paternal','CNL_LOH',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),str(p_parent_cnv)+'->'+str(p_cnv)])
                            clone.maternal_cnvs[i] = m_cnv
                            clone.paternal_cnvs[i] = p_cnv
                            clone.changes[i] = 'CNL_LOH'
                            continue
                        
                        # check mirrored CNV
                        if clone.changes and clone.changes[i-1] in ['REGULAR', 'NONE']:
                            if m_cnv != p_cnv and m_cnv == clone.paternal_cnvs[-1] and clone.maternal_cnvs[-1] == p_cnv:
                                clone.maternal_cnvs[i] = m_cnv
                                clone.paternal_cnvs[i] = p_cnv
                                if clone.changes[i-1] == 'REGULAR':
                                    changes.pop()
                                    changes.pop()
                                changes.append([clone.parent.name,clone.name,'maternal','mirrored cnv',ref['Chromosome'][i-1]+':'+str(ref['Start'][i-1])+'-'+str(ref['End'][i-1]),str(clone.parent.maternal_cnvs[i-1])+'->'+str(clone.maternal_cnvs[i-1])])
                                changes.append([clone.parent.name,clone.name,'maternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),str(m_parent_cnv)+'->'+str(m_cnv)])
                                changes.append([clone.parent.name,clone.name,'paternal','mirrored cnv',ref['Chromosome'][i-1]+':'+str(ref['Start'][i-1])+'-'+str(ref['End'][i-1]),str(clone.parent.paternal_cnvs[i-1])+'->'+str(clone.paternal_cnvs[i-1])])
                                changes.append([clone.parent.name,clone.name,'paternal','mirrored cnv',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),str(p_parent_cnv)+'->'+str(p_cnv)])
                                clone.changes[i-1] = 'mirrored cnv'
                                clone.changes[i] = 'mirrored cnv'
                                continue

                        if m_cnv != m_parent_cnv:
                            if m_cnv == 0:
                                changes.append([clone.parent.name,clone.name,'maternal','del',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),str(m_parent_cnv)+'->'+str(m_cnv)])
                            else:
                                changes.append([clone.parent.name,clone.name,'maternal','dup',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),str(m_parent_cnv)+'->'+str(m_cnv)])
                        
                        if p_cnv != p_parent_cnv:
                            if p_cnv == 0:
                                changes.append([clone.parent.name,clone.name,'paternal','del',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),str(p_parent_cnv)+'->'+str(p_cnv)])
                            else:
                                changes.append([clone.parent.name,clone.name,'paternal','dup',ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),str(p_parent_cnv)+'->'+str(p_cnv)])
                        clone.maternal_cnvs[i] = m_cnv
                        clone.paternal_cnvs[i] = p_cnv
                        clone.changes[i] = 'REGULAR'
            
            ref[clone.name+'_maternal_cnvs'] = clone.maternal_cnvs
            ref[clone.name+'_paternal_cnvs'] = clone.paternal_cnvs
            queue.extend(clone.children)

        # ref.to_csv('original_cnv.csv')

        # merge segments
        merged_rows = []

        # Iterate through the DataFrame and merge rows
        current_row = None

        for index, row in ref.iterrows():
            if current_row is None:
                current_row = row.copy()
            else:
                if list(current_row.iloc[3:]) == list(row.iloc[3:]):
                    current_row['End'] = max(current_row['End'], row['End'])
                else:
                    merged_rows.append(current_row)
                    current_row = row.copy()

        # Add the last row to the merged rows
        if current_row is not None:
            merged_rows.append(current_row)

        # Create a new DataFrame with merged rows
        # merged_ref = pd.DataFrame(merged_rows).reset_index(drop=True)
        merged_ref = pd.DataFrame(merged_rows)

        return merged_ref, changes, maternal_genome, paternal_genome

    def _generate_fasta_for_each_clone(self, root, ref, changes, maternal_genome, paternal_genome, outdir):
        queue = deque(root.children)

        while queue:
            clone = queue.popleft()
            clone.maternal_fasta = os.path.join(outdir, clone.name+'_maternal.fasta')
            clone.paternal_fasta = os.path.join(outdir, clone.name+'_paternal.fasta')
            
            with open(clone.maternal_fasta, 'w') as m_output:
                with open(clone.paternal_fasta, 'w') as p_output:
                    for chrom in maternal_genome.keys():
                        m_output.write('>'+chrom+'\n')
                        p_output.write('>'+chrom+'\n')
                        chrom_ref = ref[ref['Chromosome'] == chrom]
                        for index, row in chrom_ref.iterrows():
                            m_cnv = int(row[clone.name+'_maternal_cnvs'])
                            p_cnv = int(row[clone.name+'_paternal_cnvs'])
                            start = int(row['Start'])
                            end = int(row['End'])

                            # handle CNN_LOH
                            segment = chrom + ':' + str(start) + '-' + str(end)
                            try:
                                cnv_type = changes[(changes['Child']==clone.name) & (changes['Segment']==segment)]['Type'][0]
                            except:
                                cnv_type = 'REGULAR'
                            m_sequence = maternal_genome[chrom][start-1:end]
                            p_sequence = paternal_genome[chrom][start-1:end]
                            if cnv_type == 'CNN_LOH':
                                if m_cnv != 0:
                                    new_m_cnv = random.randint(1, m_cnv -1)
                                    new_p_cnv = m_cnv - new_m_cnv
                                    cnv_m_sequence = m_sequence * new_m_cnv
                                    cnv_p_sequence = m_sequence * new_p_cnv
                                else:
                                    new_m_cnv = random.randint(1, p_cnv -1)
                                    new_p_cnv = m_cnv - new_m_cnv
                                    cnv_m_sequence = p_sequence * new_m_cnv
                                    cnv_p_sequence = p_sequence * new_p_cnv
                            elif cnv_type == 'GOH':
                                seq_len = len(m_sequence)
                                random_snp_no = random.randint(5, 30)
                                random_snps = {snp : random.sample(['A','T','C','G'], 2) for snp in random.sample(range(seq_len), random_snp_no)}
                                new_m_sequence = ''
                                new_p_sequence = ''
                                snp_pos = random_snps.keys()
                                for pos in range(len(m_sequence)):
                                    if pos in snp_pos:
                                        new_m_sequence += random_snps[pos][0]
                                        new_p_sequence += random_snps[pos][1]
                                    else:
                                        new_m_sequence += m_sequence[pos]
                                        new_p_sequence += p_sequence[pos]
                                cnv_m_sequence = new_m_sequence * m_cnv
                                cnv_p_sequence = new_p_sequence * p_cnv
                            else:   
                                cnv_m_sequence = m_sequence * m_cnv
                                cnv_p_sequence = p_sequence * p_cnv
                            m_output.write(cnv_m_sequence)
                            p_output.write(cnv_p_sequence)
                            clone.maternal_fasta_length += len(cnv_m_sequence)
                            clone.paternal_fasta_length += len(cnv_p_sequence)
                        m_output.write('\n')
                        p_output.write('\n')
            queue.extend(clone.children)

    def _out_cnv_profile(self, root, ref, changes, outdir):
        # out cnv profile csv
        df = ref[['Chromosome', 'Start', 'End']]
        queue = deque([root])
        while queue:
            clone = queue.popleft()
            df[clone.name] = ref[clone.name+'_maternal_cnvs'].astype(str) + '|' + ref[clone.name+'_paternal_cnvs'].astype(str)
            queue.extend(clone.children)
        df.to_csv(os.path.join(outdir, 'cnv_profile.csv'), index=False)

        # out maternal cnv matrix
        indexes = ref['Chromosome'] + ':' + ref['Start'].astype(str) + '-' + ref['End'].astype(str)
        m_cnv = ref.filter(like='maternal_cnvs')
        m_cnv.index = indexes
        m_cnv.to_csv(os.path.join(outdir, 'maternal_cnv_matrix.csv'))

        # out paternal cnv matrix
        p_cnv = ref.filter(like='paternal_cnvs')
        p_cnv.index = indexes
        p_cnv.to_csv(os.path.join(outdir, 'paternal_cnv_matrix.csv'))

        # out changes profile
        columns = ['Parent', 'Child', 'Haplotype', 'Type', 'Segment', 'Change']
        change_df = pd.DataFrame(data=changes, columns=columns)
        change_df.to_csv(os.path.join(outdir, 'changes.csv'), index=False)
        return change_df

    def _merge_fasta_for_each_clone(self, root, outdir):
        # merge fasta for each clone
        queue = deque([root])
        while queue:
            clone = queue.popleft()
            clone.fasta = os.path.join(outdir, clone.name+'.fasta')
            command = """sed '/^>chr/ s/$/-A/' {0} > {1} && sed '/^>chr/ s/$/-B/' {2} >> {1}""".format(clone.maternal_fasta, clone.fasta, clone.paternal_fasta)
            code = os.system(command)
            queue.extend(clone.children)

    def _wgsim_process(self, pe_reads, fasta, fq1, fq2):
        #where yyy is the read length, zzz is the error rate and $xxx * $yyy = 10000000.
        command = self.wgsim_path + " -e {0} -d {1} -s 35 -N {2} -1 {3} -2 {3} -r0 -R0 -X0 {4} {5} {6}".format(self.error_rate,self.insertion_size,pe_reads,self.reads_len,fasta,fq1,fq2)
        # logging.info(command)
        code = os.system(command)

    def _generate_fastq(self, root, outdir):
        pool = Pool(self.wgsim_thread)

        queue = deque([root])
        while queue:
            clone = queue.popleft()
            pe_reads = round((clone.maternal_fasta_length + clone.paternal_fasta_length)*self.clone_coverage/(self.reads_len*2), 6)
            fq1 = os.path.join(outdir, clone.name+'_r1.fq')
            fq2 = os.path.join(outdir, clone.name+'_r2.fq')
            clone.fq1 = fq1
            clone.fq2 = fq2

            pool.apply_async(self._wgsim_process, (pe_reads,clone.fasta,fq1,fq2))
            queue.extend(clone.children)
            
        pool.close()
        pool.join()

    # def _downsampling_fastq(self, root, outdir):
    #     queue = deque([root])
    #     assigns = utils.assign_cells_to_clones(self.cell_no, self.clone_no)
    #     while queue:
    #         clone = queue.pop()
    #         clone_dir = os.path.join(outdir, clone.name)
    #         if not os.path.exists(clone_dir):
    #             os.makedirs(clone_dir)
    #         clone.cell_no = assigns.pop()
    #         clone.ratio = round(clone.cell_no/self.cell_no, 2)
    #         total_reads_no = int(int(subprocess.check_output(['wc', '-l', clone.fq1]).split()[0])/4)
    #         per_cell_reads_no = int(total_reads_no*self.cell_coverage/self.clone_coverage)
    #         sampling_results = utils.random_sampling(total_reads_no, clone.cell_no, per_cell_reads_no)

    #         for i in range(0, clone.cell_no):
    #             lines = sampling_results[i]

    #             # generate index file
    #             index_file = os.path.join(clone_dir, clone.name + '_cell'+str(i)+'.index')
    #             with open(index_file, 'w') as output:
    #                 for line in lines:
    #                     output.write(str(line)+'\n')
    #             if self.mode ==0: # normal mode
    #                 cell_fq1 = os.path.join(clone_dir, clone.name + '_cell'+str(i)+'_r1.fq')
    #                 cell_fq2 = os.path.join(clone_dir, clone.name + '_cell'+str(i)+'_r2.fq')
    #                 command = "awk 'NR == FNR{ ind[$1]; next }(FNR in ind)' "+index_file+" "+clone.fq1+" > " + cell_fq1
    #                 code = os.system(command)
    #                 command = "awk 'NR == FNR{ ind[$1]; next }(FNR in ind)' "+index_file+" "+clone.fq2+" > " + cell_fq2
    #                 code = os.system(command)

    #                 # delete index file
    #                 os.remove(index_file)
    #             else: # 10X barcode mode
    #                 pass   
    #         os.remove(clone.fq1)
    #         os.remove(clone.fq2)                
    #         queue.extend(clone.children)

    def alignment(self):
        fastq_dir = os.path.join(self.outdir, 'fastq')
        bam_dir = os.path.join(self.outdir, 'bam')
        picard_tmp_dir = os.path.join(self.outdir, 'tmp')

        if not os.path.exists(bam_dir):
            os.makedirs(bam_dir)
        if not os.path.exists(picard_tmp_dir):
            os.makedirs(picard_tmp_dir)

        files = glob(fastq_dir+"/*.fq")
        clones = utils.get_all_clones(files)
        for clone in clones:
            fq1 = os.path.join(fastq_dir, clone + "_r1.fq")
            fq2 = os.path.join(fastq_dir, clone + "_r2.fq")
            sam_file = os.path.join(bam_dir, clone+".sam")
            bam_file = os.path.join(bam_dir, clone+".bam")
            sorted_bam_file = os.path.join(bam_dir, clone+".sorted.bam")
            dedup_bam_file = os.path.join(bam_dir, clone+".sorted.dedup.bam")
            dedup_metrics_file = os.path.join(bam_dir, clone+".sorted.dedup.metrics.txt")

            #run bwa
            logging.info('BWA alignment for {0}...'.format(clone))
            command = "{0} mem -M -t {1} {2} {3} {4} > {5}".format(self.bwa_path, self.bwa_thread, self.ref_genome, fq1, fq2, sam_file)
            code = os.system(command)

            # samtools sam to bam
            logging.info('Samtools sam to bam for {0}...'.format(clone))
            command = "{0} view -bS {1} > {2}".format(self.samtools_path, sam_file, bam_file)
            code = os.system(command)

            # # run picard sort
            # logging.info('Picard MarkDuplicates for {0}...'.format(clone))
            # command = """java -Xmx40G -Djava.io.tmpdir={3} -jar {0} SortSam \
            #             INPUT={1} OUTPUT={2} \
            #             SORT_ORDER=coordinate TMP_DIR={3}""".format(self.picard_path, bam_file, sorted_bam_file, picard_tmp_dir)
            # code = os.system(command)

            # #run samtools build index
            # logging.info('Samtools build index for {0}...'.format(clone))
            # command = "{0} index {1}".format(self.samtools_path, sorted_bam_file)
            # code = os.system(command)

            # # run picard dedup
            # logging.info('Picard MarkDuplicates for {0}...'.format(clone))
            # command = """java -Xmx40G -Djava.io.tmpdir={4} -jar {0} MarkDuplicates \
            #             REMOVE_DUPLICATES=true \
            #             I={1} O={2} \
            #             METRICS_FILE={3} \
            #             PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_VERSION=null \
            #             PROGRAM_GROUP_NAME=MarkDuplicates TMP_DIR={4}""".format(self.picard_path, sorted_bam_file, dedup_bam_file, dedup_metrics_file, picard_tmp_dir)
            # code = os.system(command)

            #  # run picard buildindex
            # logging.info('Picard BuildBamIndex for {0}...'.format(clone))
            # command = "java -Djava.io.tmpdir={4} -jar {0} BuildBamIndex I={1} -Djava.io.tmpdir={2} TMP_DIR={2}".format(self.picard_path, dedup_bam_file, picard_tmp_dir)
            # code = os.system(command)
            # break

            # clear fastq and sam file
            # os.remove(fq1)
            # os.remove(fq2)
            # os.remove(sam_file)

    def downsampling_bam(self):
        bam_dir = os.path.join(self.outdir, 'bam')
        barcodes = []
        files = glob(bam_dir+"/*.bam")
        clones = [Path(file).stem for file in files]
        assign_cells = utils.assign_cells_to_clones(self.cell_no, self.clone_no)
        cell_ratio = round(self.cell_coverage/self.clone_coverage, 2)

        for index, clone in enumerate(clones):
            clone_bam_file = os.path.join(bam_dir, clone+".bam")
            clone_cell_no = assign_cells[index]
            clone_cell_bam_dir = os.path.join(self.outdir, 'bam', clone)

            if not os.path.exists(clone_cell_bam_dir):
                os.makedirs(clone_cell_bam_dir)
            
            # run samtools to subsampling
            logging.info('Samtools downsampling bam for {0}...'.format(clone))
            for i in range(clone_cell_no):
                cell_name = clone + '_cell' + str(i+1)
                barcodes.append(cell_name)
                cell_bam_file = os.path.join(clone_cell_bam_dir, cell_name+'.bam')
                ratio = cell_ratio + i
                command = "{0} view -b -s {1} {2} > {3}".format(self.samtools_path, ratio, clone_bam_file, cell_bam_file)
                print(command)
                code = os.system(command)

        # write barcodes file
        profile_dir = os.path.join(self.outdir, 'profile')
        barcodes_file = os.path.join(profile_dir, 'barcodes.txt')
        with open(barcodes_file, 'w') as output:
            for barcode in barcodes:
                output.write(barcode+'\n')

    def process_cell_bam(self):
        bam_dir = os.path.join(self.outdir, 'bam')
        barcodes = []
        picard_tmp_dir = os.path.join(self.outdir, 'tmp')

        if not os.path.exists(picard_tmp_dir):
            os.makedirs(picard_tmp_dir)

        # write barcodes file
        profile_dir = os.path.join(self.outdir, 'profile')
        barcodes_file = os.path.join(profile_dir, 'barcodes.txt')
        with open(barcodes_file, 'r') as output:
            for line in output.readlines():
                barcodes.append(line.strip())
        
        for barcode in barcodes:
            clone = barcode.split('_')[0]
            bam_file = os.path.join(bam_dir, clone, barcode+'.bam')
            sorted_bam_file = os.path.join(bam_dir, clone, barcode+".sorted.bam")
            dedup_bam_file = os.path.join(bam_dir, clone, barcode+".sorted.dedup.bam")
            dedup_metrics_file = os.path.join(bam_dir, clone, barcode+".sorted.dedup.metrics.txt")

            # run picard sort
            logging.info('Picard MarkDuplicates for {0}...'.format(barcode))
            command = """java -Xmx40G -Djava.io.tmpdir={3} -jar {0} SortSam \
                        INPUT={1} OUTPUT={2} \
                        SORT_ORDER=coordinate TMP_DIR={3}""".format(self.picard_path, bam_file, sorted_bam_file, picard_tmp_dir)
            code = os.system(command)

            #run samtools build index
            logging.info('Samtools build index for {0}...'.format(barcode))
            command = "{0} index {1}".format(self.samtools_path, sorted_bam_file)
            code = os.system(command)

            # run picard dedup
            logging.info('Picard MarkDuplicates for {0}...'.format(barcode))
            command = """java -Xmx40G -Djava.io.tmpdir={4} -jar {0} MarkDuplicates \
                        REMOVE_DUPLICATES=true \
                        I={1} O={2} \
                        METRICS_FILE={3} \
                        PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_VERSION=null \
                        PROGRAM_GROUP_NAME=MarkDuplicates TMP_DIR={4}""".format(self.picard_path, sorted_bam_file, dedup_bam_file, dedup_metrics_file, picard_tmp_dir)
            code = os.system(command)

             # run picard buildindex
            logging.info('Picard BuildBamIndex for {0}...'.format(barcode))
            command = "java -Djava.io.tmpdir={4} -jar {0} BuildBamIndex I={1} -Djava.io.tmpdir={2} TMP_DIR={2} VALIDATION_STRINGENCY=LENIENT".format(self.picard_path, dedup_bam_file, picard_tmp_dir)
            code = os.system(command)


    def sim_dataset(self):
        logging.info("Start simulation process...")
        fasta_dir = os.path.join(self.outdir, 'fasta')
        fastq_dir = os.path.join(self.outdir, 'fastq')
        profile_dir = os.path.join(self.outdir, 'profile')

        # create if dir not exist
        if not os.path.exists(fasta_dir):
            os.makedirs(fasta_dir)
        if not os.path.exists(fastq_dir):
            os.makedirs(fastq_dir)
        if not os.path.exists(profile_dir):
            os.makedirs(profile_dir)

        m_fasta = os.path.join(fasta_dir, 'normal_maternal.fasta')
        p_fasta = os.path.join(fasta_dir, 'normal_paternal.fasta')
        phase_file = os.path.join(profile_dir, 'phases.csv')

        # generate normal fasta with snps 
        logging.info("Building normal fasta file...")
        self._get_chrom_sizes()
        self._buildGenome(m_fasta, p_fasta, phase_file)
        
        # generate random clone tree and set root as normal clone
        logging.info('Generating random evolution tree...')
        root = random_tree.generate_random_tree_balance(self.clone_no, self.max_cnv_tree_depth)
        root.maternal_fasta = m_fasta
        root.paternal_fasta = p_fasta
        normal_fasta_length = sum(self.chrom_sizes.values())
        root.maternal_fasta_length = normal_fasta_length
        root.paternal_fasta_length = normal_fasta_length

        # output tree graph
        logging.info('Drawing tree graph...')
        random_tree.draw_tree_to_pdf(root, os.path.join(profile_dir, 'tree.pdf'))
        

        # generate cnv for each clone
        logging.info('Generating CNV profile for each clone...')
        ref = self._split_chr_to_bins('all')
        new_ref, changes, maternal_genome, paternal_genome = self._generate_cnv_profile_for_each_clone(root, ref, m_fasta, p_fasta)
        new_changes = self._out_cnv_profile(root, new_ref, changes, profile_dir)

        # generate fasta file for each clone
        logging.info('Generating fasta file for each clone...')
        self._generate_fasta_for_each_clone(root, new_ref, new_changes, maternal_genome, paternal_genome, fasta_dir)

        # add normal to the tree
        # normal = random_tree.TreeNode('normal')
        # normal.children.append(root)
        # root.parent = normal
        # normal.maternal_fasta = m_fasta
        # normal.paternal_fasta = p_fasta
        # normal_fasta_length = sum(chrom_sizes.values())
        # normal.maternal_fasta_length = normal_fasta_length
        # normal.paternal_fasta_length = normal_fasta_length

        # merge maternal and paternal to one fasta file
        logging.info('Merging maternal and paternal fasta file for each clone...')
        self._merge_fasta_for_each_clone(root, fasta_dir)

        # generate fastq for each clone
        logging.info('Generating fastq file for each clone...')
        self._generate_fastq(root, fastq_dir)

        # sampling
        # logging.info('Generating fastq file for cells of each clone...')
        # self._downsampling_fastq(root, fastq_dir)

        # output tree newick
        tree_newick = os.path.join(profile_dir, 'tree.newick')
        result = random_tree.tree_to_newick(root)
        with open(tree_newick, 'w') as output:
            output.write(result)