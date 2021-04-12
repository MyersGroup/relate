#include "catch_amalgamated.hpp"

TEST_CASE( "Testing Mutation Cathegories "){

  //define mutation types
  std::string alphabet = "ACGT";
  std::map<std::string, int> dict_mutation_pattern;

  //A-C, A-G, A-T, C-A, C-G, C-T
  //T-G, T-C, T-A, G-T, G-C, G-A
  //plus 16 flanks per mutation type, i.e. 6*16 = 96

  std::string pattern;
  int index = 0;
  for(std::string::iterator it_str1 = alphabet.begin(); it_str1 != alphabet.end(); it_str1++){
    for(std::string::iterator it_str2 = alphabet.begin(); it_str2 != alphabet.end(); it_str2++){
      pattern = *it_str1;
      pattern += *it_str2;
      dict_mutation_pattern[pattern + "A-C"] = index;
      dict_mutation_pattern[pattern + "T-G"] = index;
      index++;
      dict_mutation_pattern[pattern + "A-G"] = index;
      dict_mutation_pattern[pattern + "T-C"] = index;
      index++;
      dict_mutation_pattern[pattern + "A-T"] = index;
      dict_mutation_pattern[pattern + "T-A"] = index;
      index++;
      dict_mutation_pattern[pattern + "C-A"] = index;
      dict_mutation_pattern[pattern + "G-T"] = index;
      index++;
      dict_mutation_pattern[pattern + "C-G"] = index;
      dict_mutation_pattern[pattern + "G-C"] = index;
      index++;
      dict_mutation_pattern[pattern + "C-T"] = index;
      dict_mutation_pattern[pattern + "G-A"] = index;
      index++;
    }   
  }
  int num_mutation_cathegories = index;
  REQUIRE(num_mutation_cathegories == 96);

}
