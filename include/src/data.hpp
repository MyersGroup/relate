#ifndef DATA_HPP
#define DATA_HPP

#include "collapsed_matrix.hpp"
#include "gzstream.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <list>
#include <vector>
#include <bitset>
#include <map>
#include <limits>
#include <tgmath.h>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <cassert>

class gzip{

  private:

    bool is_gzipped = false;
    bool is_open    = false;

  public:

   gzip(){};

   //function that takes in a file pointer, filename and opens the file
   FILE* open(const char* filename, const char* mode);
   
   //function that takes in a file pointer and closes the file.
   void close(FILE* fp);

};

//struct recording all the data needed for building anc
struct Data{

  std::string name;
  int N, L; //number of sequences, number of SNPs
  int Ne; //effective population size
  double mu; //mutation rate
  double theta, ntheta; //mutation probability for painting. set to 0.001

  CollapsedMatrix<char> sequence; //sequence matrix, containing 0 and 1
  std::vector<int> state; //vector specifying state of SNP (e.g., whether to use or not for bl estimation)
	std::vector<int> bp_pos;
	std::vector<int> dist;    //vector specifying location of each SNP along the genome
  std::vector<double> r;   //vector of recombination distances from one SNP to the next
  std::vector<double> rpos; //vector of cumulative recombination distances

  ///////////

  Data(){}
  //Constructor, which only assigns values to N, L, Ne, mu
  Data(int N, int L, int Ne = 3e4, double mu = 1.25e-8);

  //Constructor, which reads files in binary format for fast io
  Data(const char* filename_sequence, const char* filename_pos, const char* filename_dist, const char* filename_rec, const char* filename_rpos, const char* filename_state, int Ne = 3e4, double mu = 1.25e-8);
  //Constructor, which reads pos, and param files in bin format
  Data(const char* filename_dist, const char* filename_param, int Ne = 3e4, double mu = 1.25e-8);
 
  ///////////

  void MakeChunks(const std::string& filename_haps, const std::string& filename_sample, const std::string& filename_map, const std::string& filename_dist, const std::string& file_out, bool use_transition, float max_memory = 5);
  //void MakeChunks2(const std::string& filename_haps, const std::string& filename_sample, const std::string& filename_map, const std::string& filename_dist);

  ///////////

  void WriteSequenceAsBin(const char* filename);
  void ReadSequenceFromBin(const char* filename);

  //to read/write vectors from a file into v
  template<typename T> void WriteVectorAsBin(std::vector<T>& v, const char* filename){
      FILE* pf = fopen(filename, "wb");

      assert(pf != NULL);
      unsigned int size = v.size();
      fwrite(&size, sizeof(unsigned int), 1, pf);
      fwrite(&v[0], sizeof(T), size, pf);

      fclose(pf);
    }
  template<typename T> void ReadVectorFromBin(std::vector<T>& v, const char* filename){
      FILE* pf = fopen(filename, "rb");

      assert(pf != NULL);
      unsigned int size;
      fread(&size, sizeof(unsigned int), 1, pf);
      v.resize(size);
      fread(&v[0], sizeof(T), size, pf);

      fclose(pf);
  }

};

struct Element{
  int pos;
  double rate;
};

class haps{

  //class to read/write bed.
  //define with file pointer or filename and never use it to write and read.

  private:

    int N, L;
    FILE* fp;
    gzip g; 
    char* line;

  public:

    char chr[1024];
    char rsid[1024];
    char ancestral[1024], alternative[1024];

    haps(const char* filename_haps, const char* filename_sample){ 

      //fp = fopen(filename_sample, "r");
      fp = g.open(filename_sample, "r");
      assert(fp);
      N = 0;
			char id1[1024], id2[1024], dummy[1024];
			assert(fscanf(fp, "%s %s %s", id1, id2, dummy) == 3); //header
			assert(fscanf(fp, "%s %s %s", id1, id2, dummy) == 3); //header
			while(fscanf(fp, "%s %s %s", id1, id2, dummy) == 3){
        if(strcmp(id1, id2) == 0){
          N += 2;
				}else{
          N++;
				}
			}
      g.close(fp);

			//fp = fopen(filename_haps, "r");
      fp = g.open(filename_haps, "r");
      assert(fp);
      L = 0;
      while(!feof(fp)){
        if(fgetc(fp) == '\n'){
          L++;
        }
      }
      g.close(fp);

      //fp = fopen(filename_haps, "r");
      fp = g.open(filename_haps, "r");
      assert(fp);
      line = (char*) malloc(2*N+10);

    }

    void ReadSNP(std::vector<char>& sequence, int& bp); //gets hap info for SNP
    void DumpSNP(std::vector<char>& sequence, int bp, FILE* fp_out); //dumps hap info for SNP
    void CloseFile(){g.close(fp);};

    int GetN(){return(N);}
    int GetL(){return(L);}

};

class map{

  private:

    FILE* fp;
    char buffer[30];

  public:
 
    std::vector<int> bp;
    std::vector<double> gen_pos;
    map(const char* filename);

};

struct fasta{

  std::string seq;
  void Read(const std::string filename);

};

#endif //DATA_HPP
