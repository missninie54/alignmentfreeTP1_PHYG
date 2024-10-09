
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

def encode_nucl(letter):
    """ Transforms a nucleotide into its numerical equivalent. 
    :param str letter: The nucleotide read.
    :return int: Numerical representation of the nucleotide.
    """
    encoding = {'A': 0, 'C' : 1, 'T' : 2, 'G' : 3}
    return encoding[letter]

def encode_kmer(seq, k):
    """ Transforms a k-mer into an integer with binary values of nucleotides.
    :param str seq: The sequence of nucleotides.
    :param int k: The length of the k-mer.
    :return int: Encoded k-mer.
    """
    kmer = 0
    for letter in seq[:k]:
        kmer <<= 2
        kmer += encode_nucl(letter)
    return kmer 

def enumerate_kmers(seq, k):
    """ Generates all k-mers of length k from the given sequence.
    :param str seq: The input sequence.
    :param int k: The length of each k-mer.
    :yield: Each k-mer in the sequence.
    """
    mask = (1 << (2 * k)) - 1
    kmer = encode_kmer(seq[:k], k)
    yield kmer

    for i in range(1, len(seq) - k + 1):
        kmer = ((kmer << 2) & mask) + encode_nucl(seq[i + k - 1])
        yield kmer

def jaccard(sequencesA, sequencesB, k):
    """ Calculate the Jaccard index between two lists of sequences.
    :param list sequencesA: List of sequences for organism A.
    :param list sequencesB: List of sequences for organism B.
    :param int k: Length of k-mers.
    :return float: Jaccard similarity index.
    """
    kmersA = set()
    kmersB = set()

    for seq in sequencesA:
        kmersA.update(enumerate_kmers(seq, k))
    
    for seq in sequencesB:
        kmersB.update(enumerate_kmers(seq, k))

    # Calculate intersection and union
    intersect = kmersA & kmersB
    union = kmersA | kmersB
    
    # Results
    print(f"Taille k-mers en commun pour les deux organismes : {len(intersect)}")
    print(f"Taille k-mers uniques pour les deux organismes : {len(union)}")
    
    # Calculate Jaccard index
    j = len(intersect) / len(union) if union else 0
    return j