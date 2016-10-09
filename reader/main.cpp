// C/C++ headers
#include <stdint.h>   // int64_t
#include <cstdio>     // sscanf()
#include <cstdlib>    // strtol
#include <ctime>      // clock(), CLOCKS_PER_SEC, clock_t
#include <exception>  // exception
#include <iomanip>    // setprecision()
#include <iostream>   // cout, endl
#include <new>        // bad_alloc
#include <string>     // string

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"
#include "mesh.hpp"
#include "parameter_input.hpp"
#include "io_wrapper.hpp"

//--------------------------------------------------------------------------------------
//!   \fn int main(int argc, char *argv[])
//    \brief Athena++ main program

int main(int argc, char *argv[])
{
  char *restart_filename=NULL;
  int res_flag=0;   // set to 1 if -r argument is on cmdline

  Globals::my_rank = 0;
  Globals::nranks  = 1;

// 1. read the restart filename from -r argument.
  int i=1;
  if(*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0'){
    switch(*(argv[i]+1)) {
    case 'r':                      // -r <restart_file>
      res_flag = 1;
      restart_filename = argv[++i];
      break;
    default:
      if(Globals::my_rank==0)
        std::cout<<"Usage: "<<argv[0]<<" -r rst_filename"<<std::endl;
      return(0);
      break;
    }
  }

  if(restart_filename==NULL) {
    // no input file is given
    std::cout << "### FATAL ERROR in main" << std::endl
              << "No input file or restart file is specified." << std::endl;
  }

// 2. Construct object to store input parameters, then parse input file and command line.
  ParameterInput *pinput;
  IOWrapper restartfile;
  try {
    pinput = new ParameterInput;
    if(res_flag==1) {
      restartfile.Open(restart_filename,IO_WRAPPER_READ_MODE);
      pinput->LoadFromFile(restartfile);
      // leave the restart file open for later use
    }
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed initializing class ParameterInput: "
              << ba.what() << std::endl;
    if(res_flag==1) restartfile.Close();
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
    if(res_flag==1) restartfile.Close();
    return(0);
  }


//--- Step 4. --------------------------------------------------------------------------
// Construct and initialize Reader

  Mesh *pmesh;
  try {
    pmesh = new Mesh(pinput, restartfile);
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed initializing class Reader: "
              << ba.what() << std::endl;
    if(res_flag==1) restartfile.Close();
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
    if(res_flag==1) restartfile.Close();
    return(0);
  }
  if(res_flag==1) restartfile.Close(); // close the restart file here


  delete pinput;
  delete pmesh;

#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif

  return(0);
}
