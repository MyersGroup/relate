#ifndef FILESYSTEM_HPP
#define FILESYSTEM_HPP

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <cassert>

#include <sys/types.h> // required for stat.h
#include <sys/stat.h> // no clue why required -- man pages say so

//code create/remove empty directories, remove files.

class filesys{

  private:

    struct stat info;

  public:

    filesys(){};

    void MakeDir(const char* sPath);
    void RmDir(const char* sPath);

};

#endif //FILESYSTEM_HPP
