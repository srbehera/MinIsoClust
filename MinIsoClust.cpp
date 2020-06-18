/* min-hash-k: Program to estimate Jaccard similarity using MinHash with k hash functions.
   James S. Plank
   CS494/CS594 - Advanced Algorithms and Programming
   October, 2017
 */
#include <stack> 
#include <omp.h>
#include "edlib.h"
#include <unordered_map>
#include "metrohash64.cpp"
#include <string>
#include <cstring>
#include <stdint.h>
#include <sys/time.h>
#include <fstream>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <openssl/md5.h>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <stdlib.h>
KSEQ_INIT(gzFile, gzread)

using namespace std;

void usage(const string s)
{
    fprintf(stderr, "usage: ./minhash <k>  <fasta> <kmerLen> <bucketsize> <min_element_cluster> <threshold_similarity>\n");
    if (s != "") fprintf(stderr, "%s\n", s.c_str()); 
    exit(1);
}

int main(int argc, char **argv)
{
  vector <string> seqnames;             // Sequence names
  vector <string> seqns;                // Sequences
  vector <unsigned char *> min_hashes;  // The minimum hashes for each file.
  int k;                                // The number of hashes
  int bbh;                              // Bytes per hash
  int hash_buf_size;                    // Size of the hash buffers (k*bbh) padded to 16
  unsigned int ff;                      // An integer that holds 0xffffffff
  unsigned char *hash;                  // Where we calculate the hashes for each string.
  ifstream f;
  string s1;
  int findex;
  int i, j, p, sz;
  double Intersection;                  
  int kmerLen;                          // q-gram size
  int bsize;                            // number of buckets

  /* Read the command line arguments. */

  if (argc <= 1) usage("");

  k = atoi(argv[1]); // size of minhash 
                 
  if (k <= 0) usage("k must be a number > 0");
  //bbh = atoi(argv[2]);
  //if (bbh <= 0) usage ("bbh must be a number > 0\n");
  // for (i = 3; i < argc; i++) files.push_back(argv[i]);
  
  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(argv[2], "r"); //file name
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    string sn(seq->name.s);
    seqnames.push_back(sn);
    string seqn(seq->seq.s);
    seqns.push_back(seqn);
  }
  kseq_destroy(seq);
  gzclose(fp);

  /* Calculate the number of bytes for all of the hashes, and allocate a
     hash buffer for temporary use, and a hash buffer for each data set
     to hold the minimum hashes for each data set.   Set each byte of 
     these buffers to 0xff, which is their maximum value, regardless of the
     size of the hash. */
  
  cout << "size: " << seqns.size() << endl;

  map<uint64_t, vector<string>> HASH_to_STRING;

  kmerLen = atoi(argv[3]); // q-gram size
  
  map<string, uint64_t> hashV; 
  vector<vector<uint64_t>> MHASH;
  
  /* Read the data sets.  For each value, you're going to calculate the k hashes
     and then update the minimum hashes for the data set. */

     string s = "";
  
  for (findex = 0; findex < seqns.size(); findex++) {
     string ss = seqns.at(findex);
     map<string, int> ugram;
     vector<uint64_t> mmhash;
     uint64_t hash1 = 0;
 
     for(int t = 0; t < (ss.length()-kmerLen+1); t++){
       s = ss.substr(t, kmerLen); //seqns.at(findex);
       if(ugram.find(s) == ugram.end()) {
         ugram.insert(pair<string, int>(s, 0)); 
         s = s+to_string(0); 
       }
       else { 
         ugram[s]++;
         s = s + to_string(ugram[s]); 
       } 
       for(int u = 0; u < k; u++){
         hash1 = 0;
         MetroHash64::Hash((uint8_t*)s.c_str(), s.length(), (uint8_t *)&hash1, u);
         if(t == 0) mmhash.push_back(hash1);
         else if(hash1 < mmhash.at(u)) mmhash[u] = hash1;  
       }
       if(mmhash.size() < k) cout << "error here\n";
     }
     MHASH.push_back(mmhash);
  }
  cout << MHASH.size() << endl;
  //exit(0);
  bsize = atoi(argv[4]);
  s = "";
  cout << "we are here: minhash stage completed...\n" ;
  //exit(0); 
  unordered_map<uint64_t, vector<int>> MM;

  for (i = 0; i < seqns.size(); i++) {
    for (int p = 0; p < k; p += bsize) {
      s = "";
      for(int x = p; x < (p+bsize); x++){
        s = s + to_string(MHASH[i].at(x));
      }
      uint64_t hash1 = 0;
      MetroHash64::Hash((uint8_t*)s.c_str(), s.length(), (uint8_t *)&hash1, 0);
      if(MM.find(hash1) != MM.end()) {
        if(find(MM[hash1].begin(), MM[hash1].end(), i) == MM[hash1].end()) 
          MM[hash1].push_back(i); 
      }
      else {
        vector<int> v; 
        v.push_back(i); 
        MM.insert(make_pair(hash1, v)); }
    }
  }

  for (auto& x: MM) {  
    std::vector<int>::iterator it1; 
    it1 = std::unique (x.second.begin(), x.second.end());
    x.second.resize( std::distance(x.second.begin(),it1) ); 
  }
 
  /* 
  int test_count = 0;
  for (auto& x: MM) {
    int l1 = seqns.at(i).length();
    vector<int> new_vec;
    for(int d=0; d<x.second.size(); d++) new_vec.push_back(x.second.at(d));
    //new_vec = x.second;
    cout << "\r" << test_count << " processing .. " << flush;
    for(int d=0; d<new_vec.size(); d++){
      EdlibAlignResult result;
      int indx = new_vec.at(d);
      int l2 = seqns.at(indx).length();
      if(l1 < l2) result = edlibAlign(seqns.at(i).c_str(), l1, seqns.at(indx).c_str(), l2, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
      else result = edlibAlign(seqns.at(indx).c_str(), l2, seqns.at(i).c_str(), l1, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0)); 
      double per_sim = ((double)(result.alignmentLength-result.editDistance)/(double)result.alignmentLength);
      edlibFreeAlignResult(result);
      if(per_sim < atof(argv[6]) { x.second.erase(x.second.begin()+d); } 
    }
    test_count++;
  }
  */
  cout << "we are here: bucketing completed...\n" ;
  //exit(0);
 
  ofstream myfile1;
  myfile1.open ("example_minMetro.txt");
  std::set< std::vector<int> > CLST1;
  cout << MM.size() << endl;
  int idd = 0;
 
  for(auto it = MM.begin(); it != MM.end(); ++it){
    /*  std::sort(it->second.begin(), it->second.end());
        auto last = std::unique(it->second.begin(), it->second.end());
        it->second.erase(last, it->second.end()); 
    */
 
    if(CLST1.find(it->second) == CLST1.end()) {
      CLST1.insert(it->second);
      myfile1 << "< ";
      
      for(int j=0; j<it->second.size(); j++)
        myfile1 << it->second[j] << " ";
      myfile1 << " >\n";

      if(it->second.size() > 500) {
        cout << "it->second.size():(" << (idd) << ")" << it->second.size() << endl;
        idd++;
      }
    }
  }
  myfile1.close();
  cout << "No. of clusters: " << CLST1.size() << endl;
  //exit(0); 

  vector<set<int>> CLST2;
  int ind = 0;
  int *arr = new int[seqns.size()]{0};
  
  cout << "No. of clusters: " << CLST1.size() << endl;
  
  myfile1.open ("similar.txt");
  ofstream myfile3;
  myfile3.open ("sim_cluster_only.txt");
  ofstream myfile4;
  myfile4.open ("sim_cluster_only_senames.txt");
  ofstream myfile5;
  myfile5.open ("final_contigs.faa");
  ofstream myfile6;
  myfile6.open ("final_contigs_residue.faa");
  int cid = 1, cid1 = 1;
  s = "";
  set< set<int> > SSS;
  for (i = 0; i < seqns.size(); i++) {
    int ccc = 0;
    if(arr[i] == 0){
      stack <int> st;
      st.push(i); 
      set <int> sseq;
      cout << "\r" << i << " processing.." << flush;
      myfile1 << i << "(" << seqns.at(i).length() <<"): < ";
      while(!st.empty()){
        int elm = st.top(); 
        st.pop(); 
        if(arr[elm] == 0) {
          if(elm == i) {sseq.insert(elm); arr[elm] = 1; ccc++;} 
          myfile1 << elm << "(" << seqns.at(elm).length() << ") "; 
          for (int p = 0; p < k; p += bsize) {
            s = ""; 
            for(int x = p; x < (p+bsize); x++){
              s = s + to_string(MHASH[elm].at(x)); 
            }
            uint64_t hash1 = 0;
            MetroHash64::Hash((uint8_t*)s.c_str(), s.length(), (uint8_t *)&hash1, 0);
            if(MM[hash1].size() <= 500){
              for(int j=0; j<MM[hash1].size(); j++){
                int jj = MM[hash1][j];
                if( arr[jj] == 0 ){
                  EdlibAlignResult result;
                  double per_sim = 0.0;
                  int l1 = seqns.at(i).length(), l2 = seqns.at(jj).length();
                  if(l1 < l2) 
                    result = edlibAlign(seqns.at(i).c_str(), l1, seqns.at(jj).c_str(), l2, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                  else
                    result = edlibAlign(seqns.at(jj).c_str(), l2, seqns.at(i).c_str(), l1, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                  per_sim = ((double)(result.alignmentLength-result.editDistance)/(double)result.alignmentLength);
                  edlibFreeAlignResult(result);
                  if(per_sim > atof(argv[6])) { st.push(jj); arr[jj] = 1; ccc++; sseq.insert(jj); }
                }
              }
            }
          }
        }
      }
      myfile1  << ">\n ";
      if(SSS.find(sseq) == SSS.end()){
        SSS.insert(sseq);
        if(sseq.size() > atoi(argv[5]) ) { 
          std::set<int>::iterator it; 
          myfile3 <<"< "; 
          myfile4 <<"< "; 
          int yyy = 0;
          for(it = sseq.begin(); it != sseq.end(); ++it){ 
            int sid = *it; 
            //cout << sid << "(" << seqns.at(sid).length() << ") "; 
            //cout << seqnames.at(sid) << "(" << seqns.at(sid).length() << ") ";
            myfile3 << sid << "(" << seqns.at(sid).length() << ") "; 
            myfile4 << (seqnames.at(sid)) << "(" << seqns.at(sid).length() << ") "; 
            myfile5 << ">" << cid << "_" << yyy << ":" << seqnames.at(sid) << "\n" << seqns.at(sid) << "\n";
            ++yyy; 
          } 
          myfile3 <<">\n"; 
          myfile4 <<">\n";
          ++cid;
          //exit(0); 
        }
        else {
          std::set<int>::iterator it;
          int yyy = 0;
          for(it = sseq.begin(); it != sseq.end(); ++it){
            int sid = *it;
            myfile6 << ">" << cid1 << "_" << yyy << ":\n" << seqns.at(sid) << "\n";
          } 
          ++cid1;
        }
        //++cid;
      }
    }
  }
  myfile1.close();
  myfile3.close();
  myfile4.close();
  myfile5.close();myfile6.close();
  exit(0);

}

