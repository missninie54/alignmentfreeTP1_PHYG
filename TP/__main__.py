from TP.loading import load_directory
from TP.kmers import kmer2str, encode_nucl, encode_kmer, enumerate_kmers, jaccard

if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())
    for i in range(len(filenames)):
        for j in range(i+1, len(filenames)):
            # Calculate Jaccard index
            j_index = jaccard(files[filenames[i]], files[filenames[j]], k)
            
            # Results
            print(f"L'indice de jaccard pour les deux organismes est le suivant : {filenames[i]} {filenames[j]}: {j_index:.4f}")