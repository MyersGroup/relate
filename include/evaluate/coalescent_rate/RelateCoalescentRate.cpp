#include "CoalescentRateForSection.cpp"
#include "SummarizeCoalescentRateForGenome.cpp"
#include "FinalizePopulationSize.cpp"
#include "ReEstimateBranchLengths.cpp"

#include "cxxopts.hpp"
#include <string>

int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options  
  cxxopts::Options options("RelateCoalescentRate");
  options.add_options()
    ("help", "Print help.")
    ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
    ("m,mutation_rate", "Mutation rate (float).", cxxopts::value<float>())
    ("mrate", "Filename of file containing avg mutation rates.", cxxopts::value<std::string>())
    ("coal", "Filename of file containing coalescence rates.", cxxopts::value<std::string>()) 
    ("dist", "Filename of file containing dist.", cxxopts::value<std::string>())
    ("i,input", "Filename of anc and mut files without file extension.", cxxopts::value<std::string>())
    ("o,output", "Filename for updated anc and mut files without file extension.", cxxopts::value<std::string>())
    ("poplabels", "Optional: Filename of file containing population labels. If ='hap', each haplotype is in its own group.", cxxopts::value<std::string>()) 
    ("years_per_gen", "Optional: Years per generation (float). Default: 28.", cxxopts::value<float>())
    ("num_bins", "Optional: Number of bins.", cxxopts::value<int>())
    ("first_chr", "Optional: Index of fist chr", cxxopts::value<int>())
    ("last_chr", "Optional: Index of last chr", cxxopts::value<int>())
    ("num_proposals", "Optional: Number of proposals between samples in SampleBranchLengths", cxxopts::value<int>())
    ("num_samples", "Optional: Number of samples in SampleBranchLengths", cxxopts::value<int>())
		("format", "Optional: Output file format when sampling branches. 0: Newick, 1: anc. Default: 0.", cxxopts::value<int>())
    ("seed", "Seed for MCMC in branch lengths estimation.", cxxopts::value<int>());
  
  options.parse(argc, argv);

  std::string mode = options["mode"].as<std::string>();
  
  if(!mode.compare("EstimatePopulationSize")){
 
    //variable population size.
    //Do this for whole chromosome
    //The Final Finalize should be a FinalizeByGroup  
    bool help = false;
    if(!options.count("input") || !options.count("output")){
      std::cout << "Not enough arguments supplied." << std::endl;
      std::cout << "Needed: input, output. Optional: first_chr, last_chr, poplabels, years_per_gen, num_bins." << std::endl;
      help = true;
    }
    if(options.count("help") || help){
      std::cout << options.help({""}) << std::endl;
      std::cout << "Estimate population size." << std::endl;
      exit(0);
    }  

    if(options.count("first_chr") && options.count("last_chr")){
      if(options["first_chr"].as<int>() < 0 || options["last_chr"].as<int>() < 0){
        std::cerr << "Do not use negative chr indices." << std::endl;
        exit(1);
      }
      for(int chr = options["first_chr"].as<int>(); chr <= options["last_chr"].as<int>(); chr++){ 
        CoalescentRateForSection(options, chr);
      }
      SummarizeCoalescentRateForGenome(options);  
    }else{
      CoalescentRateForSection(options);
    }    

    if(options.count("poplabels")){
      if(options["poplabels"].as<std::string>() == "hap"){
        FinalizePopulationSizeByHaplotype(options);
      }else{
        FinalizePopulationSizeByGroup(options);
      }
    }else{
      FinalizePopulationSize(options);
    }

  }else if(!mode.compare("EstimateDirPopulationSize")){
 
    //variable population size.
    //Do this for whole chromosome
    //The Final Finalize should be a FinalizeByGroup  
    bool help = false;
    if(!options.count("input") || !options.count("output") || !options.count("poplabels")){
      std::cout << "Not enough arguments supplied." << std::endl;
      std::cout << "Needed: poplabels, input, output. Optional: first_chr, last_chr, years_per_gen, num_bins." << std::endl;
      help = true;
    }
    if(options.count("help") || help){
      std::cout << options.help({""}) << std::endl;
      std::cout << "Estimate population size and directional migration." << std::endl;
      exit(0);
    }  

    if(options.count("first_chr") && options.count("last_chr")){
      if(options["first_chr"].as<int>() < 0 || options["last_chr"].as<int>() < 0){
        std::cerr << "Do not use negative chr indices." << std::endl;
        exit(1);
      }
      for(int chr = options["first_chr"].as<int>(); chr <= options["last_chr"].as<int>(); chr++){ 
        CoalescentRateDir(options, chr);
      }
      SummarizeCoalescentRateForGenome(options);  
    }else{
      CoalescentRateDir(options);
    }    

    FinalizePopulationSizeDir(options);

  }else if(!mode.compare("CoalescentRateForSection")){
  
    CoalescentRateForSection(options);
  
  }else if(!mode.compare("ConditionalCoalescentRateForSection")){
 
    ConditionalCoalescenceRateForSection(options);
  
  }else if(!mode.compare("SummarizeCoalescentRateForGenome")){

    SummarizeCoalescentRateForGenome(options);

  }else if(!mode.compare("FinalizePopulationSize")){

    bool help = false;
    if(!options.count("output")){
      std::cout << "Not enough arguments supplied." << std::endl;
      std::cout << "Needed:output. Optional: poplabels." << std::endl;
      help = true;
    }
    if(options.count("help") || help){
      std::cout << options.help({""}) << std::endl;
      std::cout << "Estimate population size." << std::endl;
      exit(0);
    }  

    if(options.count("poplabels")){
      if(options["poplabels"].as<std::string>() == "hap"){
        FinalizePopulationSizeByHaplotype(options);
      }else{
        FinalizePopulationSizeByGroup(options);
      }
    }else{
      FinalizePopulationSize(options);
    }

  }else if(!mode.compare("FinalizeCoalescenceCount")){

    bool help = false;
    if(!options.count("output")){
      std::cout << "Not enough arguments supplied." << std::endl;
      std::cout << "Needed:input, output." << std::endl;
      help = true;
    }
    if(options.count("help") || help){
      std::cout << options.help({""}) << std::endl;
      std::cout << "Count number of coalescences in epoch." << std::endl;
      exit(0);
    }  

    FinalizeCoalescenceCount(options);

  }else if(!mode.compare("ReEstimateBranchLengths")){
 
    //variable population size.
    //Do this for whole chromosome
    //The Final Finalize should be a FinalizeByGroup 
   
    bool help = false;
    if(!options.count("mutation_rate") || !options.count("mrate") || !options.count("coal") || !options.count("input") || !options.count("output")){
      std::cout << "Not enough arguments supplied." << std::endl;
      std::cout << "Needed: mutation_rate, mrate, coal, input, output. Optional: dist, seed." << std::endl;
      help = true;
    }
    if(options.count("help") || help){
      std::cout << options.help({""}) << std::endl;
      std::cout << "Estimate population size." << std::endl;
      exit(0);
    }  

    ReEstimateBranchLengths(options);

  }else if(!mode.compare("SampleBranchLengths")){
 
    //variable population size.
    //Do this for whole chromosome
    //The Final Finalize should be a FinalizeByGroup 
   
    bool help = false;
    if(!options.count("mutation_rate") || !options.count("coal") || !options.count("num_samples") || !options.count("input") || !options.count("output")){
      std::cout << "Not enough arguments supplied." << std::endl;
      std::cout << "Needed: mutation_rate, coal, num_samples, input, output. Optional: dist, num_proposals, seed." << std::endl;
      help = true;
    }
    if(options.count("help") || help){
      std::cout << options.help({""}) << std::endl;
      std::cout << "Estimate population size." << std::endl;
      exit(0);
    }  

    SampleBranchLengths(options);

  }else if(!mode.compare("SampleBranchLengthsBinary")){
 
    //variable population size.
    //Do this for whole chromosome
    //The Final Finalize should be a FinalizeByGroup 
   
    bool help = false;
    if(!options.count("mutation_rate") || !options.count("coal") || !options.count("num_samples") || !options.count("input") || !options.count("output")){
      std::cout << "Not enough arguments supplied." << std::endl;
      std::cout << "Needed: mutation_rate, coal, num_samples, input, output. Optional: dist, num_proposals, seed." << std::endl;
      help = true;
    }
    if(options.count("help") || help){
      std::cout << options.help({""}) << std::endl;
      std::cout << "Estimate population size." << std::endl;
      exit(0);
    }  

    SampleBranchLengthsBinary(options);

  }else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "EstimatePopulationSize, EstimateDirPopulationSize, ReEstimateBranchLengths, CoalescentRateForSection, ConditionalCoalescentRateForSection, SummarizeCoalescentRateForGenome, FinalizePopulationSize, SampleBranchLengths, SampleBranchLengthsBinary." << std::endl;
  
  }

  bool help = false;
  if(!options.count("mode")){
    std::cout << "Not enough arguments supplied." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    exit(0);
  }

  return 0;

}
