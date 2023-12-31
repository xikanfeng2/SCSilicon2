import os
import copy
import random
import pandas as pd
import subprocess
import logging
from . import utils
from . import random_tree


logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)
pd.options.mode.chained_assignment = None

class SCSilicon2:
    def __init__(self, ref_genome, snp_file, outdir='./', clone_no=1, cell_no=2, max_cnv_tree_depth=4, bin_len=500000, HEHOratio=0.5, cnv_prob_cutoff=0.8, clone_coverage=30, cell_coverage=0.5, reads_len=150, insertion_size=350, error_rate=0.02):
        self.ref_genome = ref_genome
        self.snp_file = snp_file
        self.outdir = outdir
        self.clone_no = clone_no
        self.cell_no = cell_no
        self.max_cnv_tree_depth = max_cnv_tree_depth
        self.bin_len = bin_len
        self.HEHOratio = HEHOratio
        self.cnv_prob_cutoff = cnv_prob_cutoff
        self.clone_coverage = clone_coverage
        self.cell_coverage = cell_coverage
        self.reads_len = reads_len
        self.insertion_size = insertion_size
        self.error_rate = error_rate
        self._check_params()
        self.samples = dict.fromkeys(['cell' + str(i+1) for i in range(self.cell_no)])
        for sample in self.samples:
            self.samples[sample] = {}


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
        utils.check_exist(snp_file=self.snp_file)
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
        utils.check_between(0,1,HEHOratio=self.HEHOratio)
        utils.check_between(0,1,cnv_prob_cutoff=self.cnv_prob_cutoff)
        utils.check_positive(clone_coverage=self.clone_coverage)
        utils.check_positive(cell_coverage=self.cell_coverage)
        utils.check_int(reads_len=self.reads_len)
        utils.check_positive(reads_len=self.reads_len)
        utils.check_int(insertion_size=self.insertion_size)
        utils.check_positive(insertion_size=self.insertion_size)
        utils.check_between(0,1,error_rate=self.error_rate)

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
        
        HEHOratio: float, optional, default: 0.5
            Ratio of heterozygous SNPs
        
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
            self.cell_no = params['ref_genome']
            del params['ref_genome']
        if 'snp_file' in params and params['snp_file'] != self.snp_file:
            self.snp_file = params['snp_file']
            del params['snp_file']
        if 'outdir' in params and params['outdir'] != self.outdir:
            self.snp_file = params['outdir']
            del params['outdir']
        if 'clone_no' in params and params['clone_no'] != self.clone_no:
            self.snp_file = params['clone_no']
            del params['clone_no']
        if 'cell_no' in params and params['cell_no'] != self.cell_no:
            self.snp_file = params['cell_no']
            del params['cell_no']
        if 'max_cnv_tree_depth' in params and params['max_cnv_tree_depth'] != self.max_cnv_tree_depth:
            self.snp_file = params['max_cnv_tree_depth']
            del params['max_cnv_tree_depth']
        if 'bin_len' in params and params['bin_len'] != self.bin_len:
            self.bin_len = params['bin_len']
            del params['bin_len']
        if 'HEHOratio' in params and params['HEHOratio'] != self.HEHOratio:
            self.HEHOratio = params['HEHOratio']
            del params['HEHOratio']
        if 'cnv_prob_cutoff' in params and params['cnv_prob_cutoff'] != self.cnv_prob_cutoff:
            self.HEHOratio = params['cnv_prob_cutoff']
            del params['cnv_prob_cutoff']
        if 'clone_coverage' in params and params['clone_coverage'] != self.clone_coverage:
            self.cluster_no = params['clone_coverage']
            del params['clone_coverage']
        if 'cell_coverage' in params and params['cell_coverage'] != self.cell_coverage:
            self.normal_frac = params['cell_coverage']
            del params['cell_coverage']
        if 'reads_len' in params and params['reads_len'] != self.reads_len:
            self.noise_frac = params['reads_len']
            del params['reads_len']
        if 'insertion_size' in params and params['insertion_size'] != self.insertion_size:
            self.noise_frac = params['insertion_size']
            del params['insertion_size']
        if 'error_rate' in params and params['error_rate'] != self.error_rate:
            self.noise_frac = params['error_rate']
            del params['error_rate']
        self._check_params()
        self.get_params()
        return self

    def get_params(self):
        print(vars(self))

    def _buildGenome(self, maternalFasta, paternalFasta, phaselist):
        allsnps = utils.parseSNPList(self.snp_file)
        phases = {}
        # m_genome = {}
        # p_genome = {}
        chrom_sizes = {}
        with open(self.ref_genome, 'r') as refinput:
            with open(maternalFasta, 'w') as out1:
                with open(paternalFasta, 'w') as out2:
                    chrom = None
                    snps = None
                    currentpos = 0 
                    for line in refinput:
                        line = line.strip()
                        if line.startswith('>'):
                            if chrom:
                                out1.write('\n')
                                out2.write('\n')
                            out1.write(line+'\n')
                            out2.write(line+'\n')
                            chrom = line.strip()[1:].split()[0]
                            # m_genome[chrom] = ''
                            # p_genome[chrom] = ''
                            chrom_sizes[chrom] = 0
                            snps = allsnps[chrom]
                            snppos = sorted(snps.keys())
                            currentsnppos = snppos.pop(0)
                            allele1 = snps[currentsnppos][0]
                            allele2 = snps[currentsnppos][1]
                        else:
                            linelen = len(line.strip())
                            if int(currentsnppos) > currentpos and int(currentsnppos) <= currentpos + linelen:
                                while int(currentsnppos) > currentpos and int(currentsnppos) <= currentpos + linelen:
                                    sindex = int(currentsnppos)-currentpos-1
                                    a = line[sindex]
                                    mline = line
                                    pline = line
                                    if random.random() < self.HEHOratio: #Heterozygous
                                        if random.random() < 0.5:
                                            a1 = allele1.lower() if a.islower() else allele1.upper()
                                            a2 = allele2.lower() if a.islower() else allele2.upper()
                                            phases[(chrom, currentsnppos)] = '0|1'
                                            mline = mline[:sindex]+a1+mline[sindex+1]
                                            pline = pline[:sindex]+a2+pline[sindex+1]
                                        else:
                                            a1 = allele2.lower() if a.islower() else allele2.upper()
                                            a2 = allele1.lower() if a.islower() else allele1.upper()
                                            phases[(chrom, currentsnppos)] = '1|0'
                                            mline = mline[:sindex]+a1+mline[sindex+1]
                                            pline = pline[:sindex]+a2+pline[sindex+1]
                                    else: #Homozygous
                                        a1 = allele1.lower() if a.islower() else allele1.upper()
                                        mline = mline[:sindex]+a1+mline[sindex+1]
                                        pline = pline[:sindex]+a1+mline[sindex+1]
                                    if snppos:
                                        currentsnppos = snppos.pop(0)
                                        allele1 = snps[currentsnppos][0]
                                        allele2 = snps[currentsnppos][1]
                                    else:
                                        break
                                # m_genome[chrom] += mline.strip()
                                # p_genome[chrom] += pline.strip()
                                chrom_sizes[chrom] += len(mline)
                                out1.write(mline)
                                out2.write(pline)
                            else:
                                # m_genome[chrom] += line.strip()
                                # p_genome[chrom] += line.strip()
                                chrom_sizes[chrom] += len(line)
                                out1.write(line)
                                out2.write(line)
                        currentpos += len(line)
        with open(phaselist, 'w') as output:
            for g in sorted(phases.keys(), key=(lambda x : (int(''.join([l for l in x[0] if l.isdigit()])), x[1]))):
                output.write('{},{},{}\n'.format(g[0], g[1], phases[g]))

        return chrom_sizes

    def _split_chr_to_bins(self, chrom_sizes, chrom):
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
            chrom_size = chrom_sizes[chrom]
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
            for chrom, chrom_size in chrom_sizes.items():
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

    def _generate_cnv_profile_for_each_clone(self, root, ref):
        stack = [root]
        changes = []
        cutoff = self.cnv_prob_cutoff
        while stack:
            clone = stack.pop()
            if clone.parent == None: # root node clone0
                for i in range(ref.shape[0]):
                    if random.random() > cutoff: # 20% cnv
                        m_cnv = utils.random_cnv()
                        p_cnv = utils.random_cnv()
                        clone.maternal_cnvs.append(m_cnv)
                        clone.paternal_cnvs.append(p_cnv)
                        if m_cnv != 1:
                            if m_cnv == 0:
                                mtype = 'del'
                            else:
                                mtype = 'dup'
                            changes.append(['normal',clone.name,'maternal',mtype,ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(m_cnv)])
                        if p_cnv != 1:
                            if p_cnv == 0:
                                ptype = 'del'
                            else:
                                ptype = 'dup'
                            changes.append(['normal',clone.name,'paternal',ptype,ref['Chromosome'][i]+':'+str(ref['Start'][i])+'-'+str(ref['End'][i]),'1->'+str(p_cnv)])
                    else: # normal
                        clone.maternal_cnvs.append(1)
                        clone.paternal_cnvs.append(1)
            else:
                # first inherit from parent
                clone.maternal_cnvs = copy.deepcopy(clone.parent.maternal_cnvs)
                clone.paternal_cnvs = copy.deepcopy(clone.parent.paternal_cnvs)

                for index, item in enumerate(clone.paternal_cnvs):
                    if random.random() > cutoff: # 20% cnv
                        m_parent_cnv = clone.maternal_cnvs[index]
                        p_parent_cnv = clone.maternal_cnvs[index]
                        if m_parent_cnv == 0:
                            m_cnv = 0
                        else:
                            m_cnv = random.randint(m_parent_cnv, 5)
                        if m_cnv != m_parent_cnv:
                            if m_cnv == 0:
                                changes.append([clone.parent.name,clone.name,'maternal','del',ref['Chromosome'][index]+':'+str(ref['Start'][index])+'-'+str(ref['End'][index]),str(m_parent_cnv)+'->'+str(m_cnv)])
                            else:
                                changes.append([clone.parent.name,clone.name,'maternal','dup',ref['Chromosome'][index]+':'+str(ref['Start'][index])+'-'+str(ref['End'][index]),str(m_parent_cnv)+'->'+str(m_cnv)])
                        if p_parent_cnv == 0:
                            p_cnv = 0
                        else:
                            p_cnv = random.randint(p_parent_cnv, 5)
                        if p_cnv != p_parent_cnv:
                            if p_cnv == 0:
                                changes.append([clone.parent.name,clone.name,'paternal','del',ref['Chromosome'][index]+':'+str(ref['Start'][index])+'-'+str(ref['End'][index]),str(p_parent_cnv)+'->'+str(p_cnv)])
                            else:
                                changes.append([clone.parent.name,clone.name,'paternal','dup',ref['Chromosome'][index]+':'+str(ref['Start'][index])+'-'+str(ref['End'][index]),str(p_parent_cnv)+'->'+str(p_cnv)])
                        clone.maternal_cnvs[index] = m_cnv
                        clone.paternal_cnvs[index] = p_cnv
            ref[clone.name+'_maternal_cnvs'] = clone.maternal_cnvs
            ref[clone.name+'_paternal_cnvs'] = clone.paternal_cnvs
            stack.extend(clone.children)

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

        return merged_ref, changes

    def _generate_fasta_for_each_clone(self, root, ref, m_fasta, p_fasta, outdir):
        stack = [root]
        maternal_genome = {}
        paternal_genome = {}

        while stack:
            clone = stack.pop()
            clone.maternal_fasta = os.path.join(outdir, clone.name+'_maternal.fasta')
            clone.paternal_fasta = os.path.join(outdir, clone.name+'_paternal.fasta')
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
                            m_sequence = maternal_genome[chrom][start-1:end] * m_cnv
                            p_sequence = paternal_genome[chrom][start-1:end] * p_cnv
                            m_output.write(m_sequence)
                            p_output.write(p_sequence)
                            clone.maternal_fasta_length += len(m_sequence)
                            clone.paternal_fasta_length += len(p_sequence)
                        m_output.write('\n')
                        p_output.write('\n')
            stack.extend(clone.children)

    def _out_cnv_profile(self, root, ref, changes, outdir):
        # out cnv profile csv
        df = ref[['Chromosome', 'Start', 'End']]
        stack = [root]
        while stack:
            clone = stack.pop()
            df[clone.name] = ref[clone.name+'_maternal_cnvs'].astype(str) + '|' + ref[clone.name+'_paternal_cnvs'].astype(str)
            stack.extend(clone.children)
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

    def _merge_fasta_for_each_clone(self, root, outdir):
        # merge fasta for each clone
        stack = [root]
        while stack:
            clone = stack.pop()
            clone.fasta = os.path.join(outdir, clone.name+'.fasta')
            command = """sed '/^>chr/ s/$/-A/' {0} > {1} && sed '/^>chr/ s/$/-B/' {2} >> {1}""".format(clone.maternal_fasta, clone.fasta, clone.paternal_fasta)
            code = os.system(command)
            stack.extend(clone.children)

    def _generate_fastq(self, root, outdir):
        stack = [root]
        while stack:
            clone = stack.pop()
            pe_reads = round((clone.maternal_fasta_length + clone.paternal_fasta_length)*self.clone_coverage/(self.reads_len*2), 6)
            fq1 = os.path.join(outdir, clone.name+'_r1.fq')
            fq2 = os.path.join(outdir, clone.name+'_r2.fq')
            clone.fq1 = fq1
            clone.fq2 = fq2
            #where yyy is the read length, zzz is the error rate and $xxx * $yyy = 10000000.
            command = "wgsim -e {0} -d {1} -s 35 -N {2} -1 {3} -2 {3} -r0 -R0 -X0 {4} {5} {6}".format(self.error_rate,self.insertion_size,pe_reads,self.reads_len,clone.fasta,fq1,fq2)
            # logging.info(command)
            code = os.system(command)
            stack.extend(clone.children)

    def _downsampling_fastq(self, root, outdir):
        stack = [root]
        assigns = utils.assign_cells_to_clones(self.cell_no, self.clone_no+1) # 1 for normal
        while stack:
            clone = stack.pop()
            clone_dir = os.path.join(outdir, clone.name)
            if not os.path.exists(clone_dir):
                os.makedirs(clone_dir)
            clone.cell_no = assigns.pop()
            clone.ratio = round(clone.cell_no/self.cell_no, 2)
            total_reads_no = int(int(subprocess.check_output(['wc', '-l', clone.fq1]).split()[0])/4)
            per_cell_reads_no = int(total_reads_no*self.cell_coverage/self.clone_coverage)
            sampling_results = utils.random_sampling(total_reads_no, clone.cell_no, per_cell_reads_no)

            for i in range(0, clone.cell_no):
                lines = sampling_results[i]

                # generate index file
                index_file = os.path.join(clone_dir, 'cell'+str(i)+'.index')
                with open(index_file, 'w') as output:
                    for line in lines:
                        output.write(str(line)+'\n')

                cell_fq1 = os.path.join(clone_dir, 'cell'+str(i)+'_r1.fq')
                cell_fq2 = os.path.join(clone_dir, 'cell'+str(i)+'_r2.fq')
                command = "awk 'NR == FNR{ ind[$1]; next }(FNR in ind)' "+index_file+" "+clone.fq1+" > " + cell_fq1
                code = os.system(command)
                command = "awk 'NR == FNR{ ind[$1]; next }(FNR in ind)' "+index_file+" "+clone.fq2+" > " + cell_fq2
                code = os.system(command)

                # delete index file
                os.remove(index_file)
                
                # lines = copy.deepcopy(sampling_results[i])
                # with open(clone.fq1, 'r') as input:
                #     with open(os.path.join(clone_dir, 'cell'+str(i)+'_r1.fq'), 'w') as output:
                #         index = 0
                #         current_line = lines.pop(0)
                #         for line in input:
                #             if index == current_line:
                #                 output.write(line)
                #                 if len(lines) == 0:
                #                     break
                #                 current_line = lines.pop(0)
                #             index += 1
                
                # lines = copy.deepcopy(sampling_results[i])
                # with open(clone.fq2, 'r') as input:
                #     with open(os.path.join(clone_dir, 'cell'+str(i)+'_r2.fq'), 'w') as output:
                #         index = 0
                #         current_line = lines.pop(0)
                #         for line in input:
                #             if index == current_line:
                #                 output.write(line)
                #                 if len(lines) == 0:
                #                     break
                #                 current_line = lines.pop(0)
                #             index += 1
                            
            stack.extend(clone.children)

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
        chrom_sizes = self._buildGenome(m_fasta, p_fasta, phase_file)
        
        # generate random clone tree
        logging.info('Generating random evolution tree...')
        root = random_tree.generate_random_tree_balance(self.clone_no, self.max_cnv_tree_depth)
        

        # generate cnv for each clone
        logging.info('Generating CNV profile for each clone...')
        ref = self._split_chr_to_bins(chrom_sizes, 'all')
        new_ref, changes = self._generate_cnv_profile_for_each_clone(root, ref)
        self._out_cnv_profile(root, new_ref, changes, profile_dir)

        # generate fasta file for each clone
        logging.info('Generating fasta file for each clone...')
        self._generate_fasta_for_each_clone(root, new_ref, m_fasta, p_fasta, fasta_dir)

        # add normal to the tree
        normal = random_tree.TreeNode('normal')
        normal.children.append(root)
        root.parent = normal
        normal.maternal_fasta = m_fasta
        normal.paternal_fasta = p_fasta
        normal_fasta_length = sum(chrom_sizes.values())
        normal.maternal_fasta_length = normal_fasta_length
        normal.paternal_fasta_length = normal_fasta_length

        # merge maternal and paternal to one fasta file
        logging.info('Merging maternal and paternal fasta file for each clone...')
        self._merge_fasta_for_each_clone(normal, fasta_dir)

        # generate fastq for each clone
        logging.info('Generating fastq file for each clone...')
        self._generate_fastq(normal, fastq_dir)

        # sampling
        logging.info('Generating fastq file for cells of each clone...')
        self._downsampling_fastq(normal, fastq_dir)

        # output tree graph and newick
        logging.info('Drawing tree graph...')
        random_tree.draw_tree_to_pdf(normal, os.path.join(profile_dir, 'tree.pdf'))
        tree_newick = os.path.join(profile_dir, 'tree.newick')
        result = random_tree.tree_to_newick(normal)
        with open(tree_newick, 'w') as output:
            output.write(result)