# Compiling and Running MinIsoClust

    1. git clone https://github.com/srbehera/MinIsoClust.git
    2. cd MinIsoClust
    3. g++ -o  MinIsoClust  MinIsoClust.cpp ./edlib/edlib/src/edlib.cpp  -std=c++11 -O3 -march=native -lz -lssl -lcrypto -fopenmp -I ./edlib/edlib/include/
    ./MinIsoClust <N> <fasta> <q> <b> <l> <th>
            N: minhash signature size 
            fasta: input fasta file name
            q: q-gram size
            b: bucket size
            l: minimum sequences in each cluster
            th: edit distance threshold
            
            
    4. sh clust_script.sh [ This command generates the final Output file]

# Output
    clst.min.tsv is the output file
