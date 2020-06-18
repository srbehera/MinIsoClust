# Running MinIsoClust

    git clone https://github.com/srbehera/MinIsoClust.git
    cd MinIsoClust
    g++ -o  MinIsoClust  MinIsoClust.cpp ./edlib/edlib/src/edlib.cpp  -std=c++11 -O3 -march=native -lz -lssl -lcrypto -fopenmp -I ./edlib/edlib/include/
    sh clust_script.sh 

# Output
    clst.min.tsv is the output file
