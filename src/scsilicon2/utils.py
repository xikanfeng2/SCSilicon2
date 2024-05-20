import numbers
import os
import random
import ntpath

def check_exist(**params):
    """Check that files are exist as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if not os.path.exists(params[p]):
            raise ValueError(
                "{} file or directory {} does not exist.".format(p, params[p]))

def check_positive(**params):
    """Check that parameters are positive as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] <= 0:
            raise ValueError(
                "Expected {} > 0, got {}".format(p, params[p]))

def check_lt_zero(**params):
    """Check that parameters are larger than zero as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] < 0:
            raise ValueError(
                "Expected {} > 0, got {}".format(p, params[p]))


def check_int(**params):
    """Check that parameters are integers as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if not isinstance(params[p], numbers.Integral):
            raise ValueError(
                "Expected {} integer, got {}".format(p, params[p]))


def check_bool(**params):
    """Check that parameters are bools as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] is not True and params[p] is not False:
            raise ValueError(
                "Expected {} boolean, got {}".format(p, params[p]))


def check_between(v_min, v_max, **params):
    """Checks parameters are in a specified range

    Parameters
    ----------

    v_min : float, minimum allowed value (inclusive)

    v_max : float, maximum allowed value (inclusive)

    params : object
        Named arguments, parameters to be checked

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] < v_min or params[p] > v_max:
            raise ValueError("Expected {} between {} and {}, "
                             "got {}".format(p, v_min, v_max, params[p]))

def check_in(choices, **params):
    """Checks parameters are in a list of allowed parameters
    Parameters
    ----------
    choices : array-like, accepted values
    params : object
        Named arguments, parameters to be checked
    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] not in choices:
            raise ValueError(
                "{} value {} not recognized. Choose from {}".format(
                    p, params[p], choices))

def randomSNPList(chrom_sizes, snp_ratio):
    snpList = {}
    for chrom, chrom_len in chrom_sizes.items():
        snpList[chrom] = {snp : random.sample(['A','T','C','G'], 2) for snp in random.sample(range(1, chrom_len+1), int(round(chrom_len * snp_ratio)))}
    return snpList

def parseSNPList(snpfile):
    snpList = {}
    with open(snpfile, 'r') as input:
        for line in input:
            info = line.strip().split('\t')
            chrom = info[0]
            pos = info[1]
            allele1 = info[2].upper()
            allele2 = info[3].upper()
            if allele1 not in ['A', 'C', 'G', 'T'] or allele2 not in ['A', 'C', 'G', 'T']:
                continue
            if chrom not in snpList:
                snpList[chrom] = {}
            snpList[chrom][pos] = (allele1, allele2)
    return snpList

def parseIgnoreList(ignorefile):
    ignorelist = []
    with open(ignorefile, 'r') as input:
        for line in input:
            if line != '':
                ignorelist.append(line.strip())
    return ignorelist

def random_cnv():
    # 定义数字和对应的概率e
    numbers = [0, 2, 3, 4, 5]
    probabilities = [0.05, 0.4, 0.3, 0.2, 0.05]

    # 使用choices函数进行随机选择，weights参数指定概率
    return random.choices(numbers, weights=probabilities)[0]

def random_mirrored_cnv():
    # 定义数字和对应的概率e
    numbers = [2, 3, 4, 5, 6, 7, 8, 9]
    probabilities = [0.05, 0.2, 0.3, 0.2, 0.1, 0.05, 0.05, 0.05]

    # 使用choices函数进行随机选择，weights参数指定概率
    return random.choices(numbers, weights=probabilities)[0]

def random_WGD():
    # 定义数字和对应的概率e
    numbers = [2, 3, 4, 5]
    probabilities = [0.5, 0.3, 0.15, 0.05]

    # 使用choices函数进行随机选择，weights参数指定概率
    return random.choices(numbers, weights=probabilities)[0]

def random_CNL():
    # 定义数字和对应的概率e
    numbers = [1, 2, 3, 4, 5]
    probabilities = [0.4, 0.3, 0.15, 0.1, 0.05]

    # 使用choices函数进行随机选择，weights参数指定概率
    return random.choices(numbers, weights=probabilities)[0]

def assign_cells_to_clones(cell_no, clone_no):
    # # Ensure at least one cell for each clone
    # clones = [1] * clone_no
    # remaining_cells = cell_no - clone_no

    # # Assign remaining cells to clones
    # for _ in range(remaining_cells):
    #     clone_index = random.randint(0, clone_no - 1)
    #     clones[clone_index] += 1

    # return clones
    # Initialize a list to hold the number of cells assigned to each clone
    clones = [0] * clone_no
    
    # Distribute the cells to the clones
    for i in range(cell_no):
        # Determine which clone to assign the current cell to
        clone_index = i % clone_no
        # Increment the cell count for the determined clone
        clones[clone_index] += 1
    
    return clones

def random_sampling(m, n, k):
    # m for total reads no, n for cell no, k for per reads no in cell
    if n * k < m:
        # Randomly sample without replacement
        sampled_lines = random.sample(range(m), n * k)
    else:
        # Randomly sample with replacement
        sampled_lines = [random.randint(0, m - 1) for _ in range(n * k)]

    # Reshape the sampled lines into a 2D array
    sampled_lines_2d = [sorted(sampled_lines[i:i+k]) for i in range(0, len(sampled_lines), k)]
    
    # add three lines to each line index
    new_sampled_lines_2d = []
    for lines in sampled_lines_2d:
        temp = []
        for line in lines:
            temp += [line*4+1, line*4+2, line*4+3, line*4+4]
        new_sampled_lines_2d.append(temp)
    return new_sampled_lines_2d

def get_all_clones(files):
    clones = []
    for file in files:
        clone = ntpath.basename(file).split('_')[0]
        if clone not in clones:
            clones.append(clone)

    return clones

def root_path():
    return os.path.dirname(os.path.abspath(__file__))

