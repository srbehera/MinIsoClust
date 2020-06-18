#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>
using namespace std;

int main(int argc, char **argv){
  std::ifstream file(argv[1]);
  std::string   line;
  ofstream myfile(argv[2]);
  while(std::getline(file, line))
  {
    char *line1 = (char *) malloc(sizeof(char)*line.length()+1);
    strcpy(line1,line.c_str());
    char *tok1 = strtok(line1, ":");
    char *tok2 = strtok(NULL, ":");
    char *tok3 = strtok(tok1, "_");
    string tok22(tok3);
    string tok222 = tok22.substr(1);
    string tok333(tok2);
    myfile << tok222 << "\t" << tok333 << endl;
  }
  file.close();
  myfile.close();

}
