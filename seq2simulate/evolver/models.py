import numpy, pylab
from scipy.linalg import expm
from codon_frequencies import CodonTable, FEqual
from types import NoneType
from random import uniform, choice, randint

class ValidationMixin(object):

    def _is_valid_dna(self, seq):
        valids = 'ATGCatgc-'
        for i in set(seq):
            if i not in valids: 
                return False
        return True, ''


class Codon(ValidationMixin, object):
    """
    Codon object containing triplet nucleotides (and possibly indels).
    
    Args:
        seq: the nucleotide sequence of the codon (default '---')
    """
    def __init__(self, seq='---'): 
        self.seq = seq

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, seq):
        if not isinstance(seq, str):
            raise ValueError('seq must be a str')
        #if len(seq) > 3:
        #    seq = seq[:3]    
        while len(seq) < 3:
            seq += '-'
        if not self._is_valid_dna(seq):
            raise ValueError('sequence is not valid')
        self._seq = seq

    def delete(self):
        self.seq = '-'*len(self.seq)
                    
    def __str__(self):
        return '<Codon {}>'.format(self._seq)
        
class Locus(object):
    """
    Locus object containing codon(s) for grouping mutations by a reference codon.
    
    Args:
        codons: a list of Codon objects or nucleotide strings (typically triplets).
    """

    def __init__(self, codons=[], loc=0):
        self.loc = loc
        self._codons = []
        for codon in codons:
            self.add_codon(codon)             
        self.history = []
   
    @property 
    def loc_aa(self):
        return self.loc+1

    @property 
    def codons(self):
        return self._codons

    @codons.setter 
    def codons(self, codons):
        if not isinstance(codons, list):
            raise ValueError('codons must be a list') 
        for codon in codons:
            if not isinstance(codon, Codon):
                raise ValueError('codons must contain Codons')
        self._codons = codons

    def add_codon(self, codon):
        if isinstance(codon, Codon):
            self._codons.append(codon)
        elif isinstance(codon, str):
            self._codons.append(Codon(seq=codon))
            
    def del_codon(self, codon_i):
        self._codons = [codon for i, codon in enumerate(self._codons) if i != codon_i]

    def __str__(self):
        return '<Locus {}>'.format(' '.join([codon.__str__() for codon in self.codons]))

    @property
    def sequence(self):
        """
        Returns compiled sequence string.
        """
        return ''.join([codon.seq for codon in self.codons])
        
    @property    
    def location_sequence(self):
        """
        Returns codon location string and sequence string.
        """
        sequence = [codon.seq for codon in self.codons]
        return [self.loc_aa]*len(sequence), sequence

    @property
    def mutations(self):
        mutations = []
        try:
            for mutation in self.history:
                if mutation[0] == 'deletion':
                    mutations.append('del{}{}'.format(self.loc_aa, convert_dna_to_aa(mutation[1])))
                elif mutation[0] == 'insertion':
                    mutations.append('ins{}{}'.format(self.loc_aa, convert_dna_to_aa(mutation[1])))
                else:
                    aa0 = convert_dna_to_aa(mutation[0])
                    aa1 = convert_dna_to_aa(mutation[1])
                    if aa0 != aa1:
                        mutations.append('{}{}{}'.format(aa0, self.loc_aa, aa1))
        except:
            pass
        return mutations

class Sequence(object):
    """
    Convenience class to collect a set of loci.
    """
    def __init__(self, seq):
        try:
            self.loci = parse_sequence_to_loci(seq) 
        except:
            self.loci = []

    @property
    def codons(self):
        return [j for i in [locus.codons for locus in self.loci] for j in i]

    @property
    def seq(self):
        return parse_loci_to_sequence_string(self.loci)
    
    @property
    def location_seq(self):
        return parse_loci_to_sequence(self.loci)

    @property
    def history(self):
        return [locus.history for locus in self.loci]

    @property
    def mutations(self):
        mutations = []
        for i in [locus.mutations for locus in self.loci]:
            if i != []:
                mutations.append(i)
        return mutations

def break_sequence_into_triplets(sequence):
    nfloor = len(sequence)//3.0
    nfloat = len(sequence)/3.0
    if nfloat > nfloor:
        nfloor += 1
    nfloor = int(nfloor)
    return [sequence[i*3:i*3+3] for i in range(nfloor)]

def convert_dna_to_aa(sequence, codon_table=None):
    if codon_table is None:
        codon_table = CodonTable()
    codon_strings = break_sequence_into_triplets(sequence)
    aas = []
    for i in codon_strings:
        if i == '---':
            aa = '-' 
        else:
            aa = getattr(codon_table, i)
        aas.append(aa)
    return ''.join([i for i in aas])
    
def parse_sequence_to_loci(sequence):
    codon_strings = break_sequence_into_triplets(sequence)
    loci = [Locus(codons=[j], loc=i) for i,j in enumerate(codon_strings)]
    return loci

def parse_loci_to_sequence(loci):
    locations = []
    sequences = []
    for locus in loci:
        location, sequence = locus.location_sequence
        locations.append(location)
        sequences.append(sequence)
    return [j for i in locations for j in i], [j for i in sequences for j in i]

def parse_loci_to_sequence_string(loci):
    sequences = []
    return ''.join([locus.sequence for locus in loci])

def goldman_Q(kappa=2.0, omega=1.0, codon_freq=None, scale_q=True, return_dict=False):
    """
    Generate Q matrix of instantaneous rates of replacement for each codon.
    Diagonals are set to values that cause row sums to equal 0. 

    Args:
        kappa: transition/transversion ratio
        omega: dN/dS
        codon_freq: codon equilibrium frequencies (found in codon_frequencies.py, dict)
        scale_q: scales Q so that the average rate of substitution at
            equilibrium equals 1. Branch lengths are thus expected number of nucleotide
            substitutions per codon. See Goldman (1994).
        return_dict: return a dictionary version of the matrix for convenience

    Returns:
        q: a 61x61 matrix
        qdict: q in dict form for convenience

    Refs:
        Goldman, N., & Yang, Z. (1994). A codon-based model of nucleotide
        substitution for protein-coding DNA sequences. Molecular Biology and Evolution,
        11(5), 725-736.
    """

    if not isinstance(kappa, (float, int)):
        raise ValueError('kappa must be a number')
    if kappa < 0:
        raise ValueError('kappa must be positive')
    if not isinstance(omega, (float, int)):
        raise ValueError('omega must be a number')
    if omega < 0:
        raise ValueError('omega must be positive')
    if not isinstance(codon_freq, (dict, NoneType)):
        raise ValueError('codon_freq must be a dictionary or None')
     
    if codon_freq is None:
        codon_freq = FEqual
    ct = CodonTable(stop_codons=False)
    codons = sorted(ct.codon_dict)
    n = len(codons)
    q = numpy.ones([n, n])
    idmatrix = numpy.identity(n)

    for i in range(n):
        codon1 = codons[i]
        for j in range(n):
            codon2 = codons[j] 
            q[i, j] = mutation_rate(codon1, codon2, omega=omega, kappa=kappa, codon_table=ct, codon_freq=codon_freq)

    # diagonal equals negative sum of row
    q[idmatrix==1] = 0.0
    q[idmatrix==1] = -q.sum(1)

    # scale Q so that average rate of substitution at equilibrium is 1
    if scale_q:
        qdiagonal = q.diagonal()
        pvalues = numpy.array([codon_freq[codon] for codon in codons])
        scaling_factor = 1.0/(-sum(qdiagonal*pvalues))
        q *= scaling_factor

    if return_dict:
        # convert matrix to codon-indexed dict
        qdict = dict(zip(codons, [dict(zip(codons, qq)) for qq in q]))
        return q, qdict

    return q

def convert_q_to_p(q, t=0.01, codon_table=None, return_dict=False):
    """

    Convert Q matrix (from goldman_Q()) to probabilities by solving matrix
    differential: P'(t) = QP(t) to yield P(t) = exp(Qt). P is a matrix of
    probabilities of a codon mutating to another over time. See Kubatko (2016).

    Args:
        q: Q matrix
        t: time or branch length over which to plot(if q was scaled, this is
        the expected number of nucleotide mutations per codon)
        codon_table: CodonTable instance
        return_dict: return a dictionary version of the matrix for convenience

    Refs: 
        Kubatko, L., Shah, P., Herbei, R., & Gilchrist, M. A. (2016). A codon model
        of nucleotide substitution with selection on synonymous codon usage. Molecular
        Phylogenetics and Evolution, 94(Pt A), 290-297. doi:10.1016/j.ympev.2015.08.026
    
    """
    if not isinstance(q, numpy.ndarray):
        raise ValueError('q must be a NumPy array')
    if not isinstance(t, (float, int)):
        raise ValueError('t must be a number')
    if not isinstance(codon_table, (CodonTable, NoneType)):
        raise ValueError('codon_table must be instance of CodonTable or None')
    if codon_table is None:
        codon_table = CodonTable(stop_codons=False)

    p = expm(q*t)

    # convert matrix to codon-indexed dict
    if return_dict:
        codons = sorted(codon_table.codon_dict)
        pdict = dict(zip(codons, [dict(zip(codons, pp)) for pp in p]))
        return p, pdict

    return p

def get_cumulative_p(p, codon_table=None, return_dict=False):
    """

    Returns cumulatively probabilities from P matrix, and codon triplets sorted
    by cumulative probability. This function is useful for sampling mutations.

    Args:
        p: P matrix (generated by convert_q_to_p())
        codon_table: CodonTable instance
        return_dict: return convenient cumulative p as dictionary with sorted codons

    """
    if not isinstance(p, numpy.ndarray):
        raise ValueError('p must be a NumPy array')
    if not isinstance(codon_table, (CodonTable, NoneType)):
        raise ValueError('codon_table must be instance of CodonTable or None')
    if codon_table is None:
        codon_table = CodonTable(stop_codons=False)
 
    codons = numpy.array(sorted(codon_table.codon_dict))
    p_cumsum = []
    p_codons = []
    for i, row in enumerate(p):
        row_sort_ind = row.argsort()
        p_cumsum.append(row[row_sort_ind].cumsum())
        p_codons.append(list(codons[row_sort_ind]))

    p_cumsum = numpy.array(p_cumsum)
    if return_dict:
        p_cumsum_dict = {codons[i]: {'p':row, 'codons':p_codons[i]} for i, row in enumerate(p_cumsum)}
        return p_cumsum, p_codons, p_cumsum_dict

    return p_cumsum, p_codons
    
    return numpy.array(p_cumsum), p_codons

 
      
def mutation_category(codon1, codon2, codon_table=None):
    """
    Args:
        codon1/codon2: three-letter strings

    Returns category of mutation one of:
        ['multisite'],
        ['synonymous', 'transversion],
        ['synonymous', 'transition],
        ['nonsynonymous', 'transversion],
        ['nonsynonymous', 'transition],
    """
    if not isinstance(codon1, str):
        raise ValueError('codon must be a three-letter str')
    if len(codon1) != 3:
        raise ValueError('codon must be a three-letter str')
    if not isinstance(codon2, str):
        raise ValueError('codon must be a three-letter str')
    if len(codon2) != 3:
        raise ValueError('codon must be a three-letter str')
    if not isinstance(codon_table, (CodonTable, NoneType)):
        raise ValueError('codon_table must be instance of CodonTable or None')

    if codon_table is None:
        codon_table = CodonTable(stop_codons=False)

    mscores = mutations(codon1, codon2)
    
    cat = None    

    if sum([nt > -1 for nt in mscores]) > 1:
        cat = ['multisite']
    else:
        nscore = max(mscores)

        nucleotide_mutation_type = None
        if nscore == 0:
            nucleotide_mutation_type = 'transition'
        if nscore == 1:
            nucleotide_mutation_type = 'transversion'

        aminoacid_mutation_type = None
        if getattr(codon_table, codon1) == getattr(codon_table, codon2):
            aminoacid_mutation_type = 'synonymous'
        else:
            aminoacid_mutation_type = 'nonsynonymous'

        cat = [aminoacid_mutation_type, nucleotide_mutation_type]

    return cat

def mutation_rate(codon1, codon2, 
    kappa=2.0, 
    omega=1.0, 
    codon_table=None, 
    codon_freq=None
    ):
    """
    Args:
        codon1/codon2: three-letter strings
        codon_table: CodonTable instance
        kappa: transition/transversion ratio
        omega: dN/dS
        codon_freq: codon equilibrium frequencies (found in codon_frequencies.py, dict)

    Returns instantaneous rate of mutation using the simplified Goldman model in Nielsen and Yang (1998).
    """
    if not isinstance(codon1, str):
        raise ValueError('codon must be a three-letter str')
    if len(codon1) != 3:
        raise ValueError('codon must be a three-letter str')
    if not isinstance(codon2, str):
        raise ValueError('codon must be a three-letter str')
    if len(codon2) != 3:
        raise ValueError('codon must be a three-letter str')
    if not isinstance(codon_table, CodonTable):
        raise ValueError('codon must be a three-letter str')
    if not isinstance(omega, (int, float)):
        raise ValueError('omega must be a number')
    if not isinstance(kappa, (int, float)):
        raise ValueError('kappa must be a number')
    if not isinstance(codon_freq, dict):
        raise ValueError('codon_freq must be an instance of dict')


    cat = mutation_category(codon1, codon2, codon_table=codon_table)
    rate = 0
    if cat == 'multisite':
        rate = 0 
    if cat == ['synonymous', 'transversion' ]:
        rate = codon_freq[codon2]
    if cat == ['synonymous', 'transition' ]:
        rate = kappa*codon_freq[codon2]
    if cat == ['nonsynonymous', 'transversion']:
        rate = omega*codon_freq[codon2]
    if cat == ['nonsynonymous', 'transition']:
        rate = omega*kappa*codon_freq[codon2]
    return rate 

def diff_index(s1, s2):
    index = []
    for i,j in zip(s1, s2):
        if i==j:
            index.append(0)
        else: 
            index.append(1)
    return index
     
def mutations(seq1, seq2):
    """
    Returns list of mutation values: 0 for transversion, 1 for transition, -1 for neither.

    Args:
        seq1, seq2: nucleotide strings
    """
    mut_dict = {
        'at': 1, 
        'ag': 0,
        'ac': 1,
        'ta': 1,
        'tg': 1,
        'tc': 0,
        'ga': 0,
        'gt': 1,
        'gc': 1,
        'ca': 1,
        'ct': 0,
        'cg': 1,
    }
    mutations = []
    for i,j in zip(seq1, seq2):
        mutations.append(mut_dict.get(i+j, -1))
    return mutations
            

def plot_p_over_time(q, t=10, codon='atg', codon_table=None, logscale=True):
    """
    Plot the probabilities of mutation of a particular codon. Requires Matplotlib.

    Args:
        q: Q matrix (generated by goldman_Q())
        t: time or branch length over which to plot(if q was scaled, this is
        the expected number of nucleotide mutations per codon)
        codon: codon triplet to plot
        logscale: log-scaled axes
         
    """

    if not isinstance(q, numpy.ndarray):
        raise ValueError('q must be a NumPy array')
    if not isinstance(codon_table, (CodonTable, NoneType)):
        raise ValueError('codon_table must be instance of CodonTable or None')
    if not isinstance(t, (float, int)):
        raise ValueError('t must be a number')
    if not isinstance(codon, str):
        raise ValueError('codon must be a string')
    if len(codon) != 3:
        raise ValueError('codon must be a string of length 3')
    
    if codon_table is None:
        codon_table = CodonTable(stop_codons=False)
    codons = sorted(codon_table.codon_dict)
    tarray = numpy.linspace(0, t, 30)
    ps = numpy.array([convert_q_to_p(q, t=t, codon_table=codon_table) for t in tarray])

    n_i = codons.index(codon)
    dim = numpy.sqrt(len(q))+1
    fig = pylab.figure(figsize=[25, 25])
    axs = [fig.add_subplot(dim, dim, i+1) for i in range(len(q))]
    p_codon = numpy.transpose([x[n_i] for x in ps])
    for n_j, ax in enumerate(axs):
        if logscale:
            ax.semilogy(tarray, p_codon[n_j], basey=10)
            ax.set_ylim([1e-5, 1])
        else:
            ax.plot(tarray, p_codon[n_j])
            ax.set_ylim([0, 1])
        ax.set_xticks([0, t])
        ax.set_yticks(ax.get_ylim())
        ax.set_title('{} ({})'.format(codons[n_j], getattr(codon_table, codons[n_j])))
    fig.suptitle('{} ({})'.format(codons[n_i], getattr(codon_table, codons[n_i])), size=30)
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.show() 

def get_mutation_from_cumulative_p(codon, p_cumsum_dict):
    """
    Makes a mutation call using a cumulative P matrix dictionary (generated by get_cumulative_p()).
    
    Args:
        codon: three-letter codon string
        p_cumsum_dict: cumulative P matrix dictionary
    Returns:
        mutated codon triplet
    """
    if not isinstance(codon, str):
        raise ValueError('codon must be a string')
    if len(codon) != 3:
        raise ValueError('codon must be a three-letter str')
    if not isinstance(p_cumsum_dict, dict):
        raise ValueError('p_cumsum_dict must be a NumPy array')
    codon_row = p_cumsum_dict[codon]
    probabilities = codon_row['p']
    codons = codon_row['codons']
    new_codon = codons[numpy.where(probabilities > uniform(0, 1))[0][0]]
    return new_codon

def make_indel(locus, index=0, ti_td=0.1, codon_freq=None):
    """
    Makes a codon-length indel at locus.
    
    Args:
        locus: and instance of Locus
        index: indel codon index
        ti_td: ratio of insertions to deletions 
        codon_feq: a codon frequency dictionary

    Returns a mutated locus.
    """
    indel_type = 'deletion'
    p_deletion = 1.0/(1+ti_td)
    if uniform(0, 1) > p_deletion:
        indel_type = 'insertion'

    if indel_type == 'deletion':
        locus.history.append([indel_type, locus.codons[index].seq, index]) 
        locus.codons[index].delete()
    elif indel_type == 'insertion':
        new_codon = choose_random_codon(codon_freq=codon_freq)
        new_codon = Codon(seq=new_codon)
        bef_aft =  choice(['before', 'after'])
        if bef_aft == 'after':
            index += 1 
        locus.codons = locus.codons[:index] + [new_codon] + locus.codons[index:]
        locus.history.append([indel_type, new_codon.seq, index]) 

def choose_random_codon(codon_freq=None):
    """
    Args:
        codon_freq: codon frequence dictionary
    """
    if codon_freq is None:
        codon_freq = FEqual
    if not isinstance(codon_freq, dict):
        raise ValueError('codon_freq must be a dictionary')

    codons = [i for i in codon_freq.keys()]
    cum_freq = numpy.array([i for i in codon_freq.values()]).cumsum()
    new_codon = codons[numpy.where(cum_freq > uniform(0, 1))[0][0]]
    return new_codon

def make_subs_in_locus(locus, q, t=0):
    """
    Mutates all codons in locus according to Q.
    """
    if not isinstance(locus, Locus):
        raise ValueError('locus must be an instance of Locus')
    if not isinstance(q, numpy.ndarray):
        raise ValueError('q must be a NumPy array')
    if not isinstance(t, (float, int)):
        raise ValueError('t must be a number')
    for i, codon in enumerate(locus.codons):
        old_seq = codon.seq
        make_sub_from_q(codon, q, t=t)
        new_seq = codon.seq
        if old_seq != new_seq:
            locus.history.append([old_seq, new_seq, i])
 
def make_sub_from_q(codon, q, t=0):
    """
    Makes a mutation call using a model Q matrix (generated by e.g. goldman_Q().
    
    Args:
        codon: a Codon instance
        q: Q matrix
        t: time or branch length over which to plot(if q was scaled, this is
        the expected number of nucleotide mutations per codon)
    Returns:
        mutated codon triplet
    """
    if not isinstance(codon, Codon):
        raise ValueError('codons must contain Codons')
    if codon.seq == '---':
        return 
    if not isinstance(q, numpy.ndarray):
        raise ValueError('q must be a NumPy array')
    if not isinstance(t, (float, int)):
        raise ValueError('t must be a number')
    p = convert_q_to_p(q, t=t)
    p_cumsum, p_codons, p_cumsum_dict = get_cumulative_p(p, return_dict=True)
    codon.seq = get_mutation_from_cumulative_p(codon.seq, p_cumsum_dict)

def plot_codon_hist(codons):
    """
    Plots a histogram of a list of items (codons in this instance).

    Args:
        codons: 1-dimensional list
    """
    if not isinstance(codons, list):
        raise ValueError('codons must be a list')

    codons = numpy.array(codons)
    set_codons = list(set(codons))
    codon_bins = numpy.array([sum(codons == i) for i in set_codons])
    sort_ind = numpy.argsort(codon_bins)[::-1]
    codon_bins = codon_bins[sort_ind]
    set_codons = [set_codons[i] for i in sort_ind] 
     
    fig, ax = pylab.subplots()
    ax.bar(range(len(codon_bins)), codon_bins, 0.8)
    ax.set_xticks(numpy.arange(len(codon_bins)))
    ax.set_xticklabels(set_codons, rotation='vertical')
    ax.set_ylabel('n')
    fig.show() 

def sample_model_mutation_probabilities(seq, q, t=0, n=100):
    if not isinstance(q, numpy.ndarray):
        raise ValueError('q must be a NumPy array')
    if not isinstance(t, (float, int)):
        raise ValueError('t must be a number')
    if not isinstance(n, (float, int)):
        raise ValueError('n must be a number')
    if t < 0:
        raise ValueError('t must be positive')
    if n <= 0:
        raise ValueError('n must be > 0')
    return [make_sub_from_q(Codon(seq=seq), q, t=t) for i in range(n)]


if __name__=='__main__':
    pass

    
     

