# Compiling and Running MinIsoClust

    git clone https://github.com/srbehera/MinIsoClust.git
    cd MinIsoClust
    g++ -o  MinIsoClust  MinIsoClust.cpp ./edlib/edlib/src/edlib.cpp  -std=c++11 -O3 -march=native -lz -lssl -lcrypto -fopenmp -I ./edlib/edlib/include/
    ./MinIsoClust <N> <fasta> <q> <b> <l> <th>
            N: minhash signature size 
            fasta: input fasta file name
            q: q-gram size
            b: bucket size
            l: minimum shared segment size
            th: edit distance threshold
            
            
    sh clust_script.sh 

# Output
    clst.min.tsv is the output file
