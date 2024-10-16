def kmer2str(val, k):
    """ Transform a kmer integer into its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved in the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

def encode_nucl(nucl):
    """ Encode a nucleotide into a 2-bit integer
    :param str nucl: The nucleotide to encode
    :return (int, int): The encoded nucleotide and its reverse complement
    """
    encoded = (ord(nucl) >> 1) & 0b11  # Extract the two bits of the ASCII code that represent the nucleotide
    rencoded = (encoded + 2) & 0b11  # Complement encoding with bit tricks. Avoid slow if statement.
    return encoded, rencoded

def xorshift(val):
    """Apply xorshift to the kmer value to compute its reverse complement"""
    val ^= (val << 13) & 0xFFFFFFFFFFFFFFFF
    val ^= (val >> 7) & 0xFFFFFFFFFFFFFFFF
    val ^= (val << 17) & 0xFFFFFFFFFFFFFFFF
    return val

def stream_kmers(seq, k):
    """Generate kmers and their reverse complements from a sequence"""
    kmer = 0
    rkmer = 0

    # Add the first k-1 nucleotides to the kmer and its reverse complement
    for i in range(k-1):
        nucl, rnucl = encode_nucl(seq[i])
        kmer |= nucl << (2 * (k-2-i))
        rkmer |= rnucl << (2 * i)

    mask = (1 << (2 * (k-1))) - 1

    # Yield the kmers
    for i in range(k-1, len(seq)):
        nucl, rnucl = encode_nucl(seq[i])
        # Remove the leftmost nucleotide from the kmer
        kmer &= mask
        # Shift the kmer to make space for the new nucleotide
        kmer <<= 2
        # Add the new nucleotide to the kmer
        kmer |= nucl
        # Make space for the new nucleotide in the reverse kmer (remove the rightmost nucleotide by side effect)
        rkmer >>= 2
        # Add the new nucleotide to the reverse kmer
        rkmer |= rnucl << (2 * (k-1))

        yield min(xorshift(kmer), xorshift(rkmer))

def filter_smallests(seq, k, sketch_size):
    """Filter and return the smallest kmers and rkmer.
    :param seq: Input sequence
    :param k: Size of kmer
    :param sketch_size: Number of smallest elements to retain
    :return: Sorted list of smallest kmers
    """
    sketch = [float('inf')] * sketch_size  # Initialize the sketch with infinity
    max_element = float('inf')  # Start with max_element as infinity

    # Iterate over each kmer generated from the sequence
    for value in stream_kmers(seq, k):
        if value < max_element:
            # Find the index of the current max element in the sketch
            idx = sketch.index(max_element)
            # Replace the max element with the new smaller value
            sketch[idx] = value
            # Update max_element to the current maximum in the sketch
            max_element = max(sketch)
    
    return sorted(sketch)