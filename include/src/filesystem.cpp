#include "filesystem.hpp"

void 
filesys::MakeDir(const char* sPath){
  
  //create directory called sPath/ if not existent
  //check if directory exists
  if(stat( sPath, &info ) != 0){
    //paint/ does not exist so create it.
    mode_t nMode = 0700; // UNIX style permissions
    int nError = 0;
#if defined(_WIN32)
    nError = _mkdir(sPath); // can be used on Windows
#else 
    nError = mkdir(sPath,nMode); // can be used on non-Windows
#endif
    if(nError != 0){
      std::cerr << "Could not create directory " << sPath << "." <<  std::endl;
      exit(1);
    }
  }

}

void
filesys::RmDir(const char* sPath){

  int nError = 0;
#if defined(_WIN32)
      nError = _rmdir(sPath); // can be used on Windows
#else 
      nError = rmdir(sPath); // can be used on non-Windows
#endif
    if(nError != 0){
      std::cerr << "Could not delete directory " << sPath << "." << std::endl;
      exit(1);
    }

}


