
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#ifdef USE_OPEMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "ALLOCATION.h"
#include "FADDEEVA.h"
#include "INV_INIT.h"
#include "IO.h"
#include "LM_FIT.h"
#include "LOG_ERROR.h"
#include "ME_SOLVER.h"
#include "MPI_INIT.h"
#include "RANDOM_NUMBER.h"
#include "RINPUT.h"
#include "STR.h"
#include "SVD.h"
#include "TIMER.h"
#include "FREE_RAM.h"

#if defined(USE_MPI) && defined(USE_OPEMP)
    #define PARALLEL_MPIOPENMP
#elif defined(USE_MPI)
    #define PARALLEL_MPI
#elif defined(USE_OPEMP)
    #define PARALLEL_OPENMP
#else
    #define PARALLEL_NONE
#endif

#define TAG_NPROF 1000
#define TAG_COORD 1001
#define TAG_DATA  1002
#define TAG_R_NPROF 2000
#define TAG_R_COORD 2001
#define TAG_R_DATA  2002

#define Input_Path "./input"

/*--------------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {

#if defined(PARALLEL_MPIOPENMP)
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if(provided < MPI_THREAD_FUNNELED){
      LOG_ERROR(ERR_LVL_ERROR, "main", \
          "Error: MPI does not provide required thread support\n");
    }
#elif defined(PARALLEL_MPI)
    MPI_Init(&argc, &argv);
#endif

#ifdef PARALLEL_MPI
    MPI_Status status;
#endif 

    STRUCT_MPI Mpi = {0};
    MPI_SETUP(&Mpi);
    Timer(&Mpi);
    const char *Filename = (argc>=2) ? argv[1] : Input_Path;

    if(!FILE_EXIST(Filename) && Mpi.rank==0){
      LOG_ERROR(ERR_LVL_ERROR, "main", "input doesn't exist! \n");
    }

    STRUCT_INPUT Input = {0};
    STRUCT_STK Stk = {0};
    STRUCT_LM LM = {0}; 
    STRUCT_PARA Para = {0};
    STRUCT_SUBSET Subset = {0};
    STRUCT_SUBSET SubsetRV = {0};

    RDINPUT(Filename, &Input, &Mpi);

    OPENMP_SETUP(&Mpi);
    rWavelength(&Input, &Stk, &Mpi);
    INIT_INV(&Input, &Stk, &LM, &Para, &Mpi, &Subset);

    int rcounts = 0, pid, src, terminal = -1, icache = 0;
    int valid = 0, actived = 1;
    double *ptr = NULL, *ptrRV = NULL;

    CACHE_INIT(&Input, &Mpi, &Stk);

    if(Mpi.rank == 0) LOG_WRITE("\n\n  --- start iversion ---  \n\n", \
        true, Input.verboselv>=2); 
    Subset.pcounts = 0;
    if(Mpi.size==1){
      if(Input.cache_prof){
        rProfilesll(&Input, &Stk);
        INVERSION_MULTI(&Input, &Stk, &Para, &LM, &Subset);

        WRITE_RESULT(&Input, &Subset, &Stk);
        
      }else{
        do{
          if(Input.cache[icache]){
            PixelMV(&Input, &Subset);
          }else{
            SubsetRV.coord[0] = Subset.coord[0];
            SubsetRV.coord[1] = Subset.coord[1];
            rProfile(&Input, &Stk, &Subset);
            SubsetRV.nProf = Subset.nProf;
            INVERSION_MULTI(&Input, &Stk, &Para, &LM, &Subset);
            WRITE_RESULT(&Input, &SubsetRV, &Stk);
          }
          
          icache++;
        }while(Subset.pcounts<Input.counts);
      }

#ifdef PARALLEL_MPI
    }else{

      if(Mpi.rank==0){
        if(Input.cache_prof){
          rProfilesll(&Input, &Stk);
        }
        ptr = Input.profbuff;

        for(pid=1, icache=0; pid<Mpi.size && Subset.pcounts<Input.counts; \
            icache++){

          if(!Input.cache_prof && Input.cache[icache]){
            PixelMV(&Input, &Subset);
            rcounts += Subset.nProf;
          }else{

            if(Subset.nProf > Input.counts-Subset.pcounts){ 
              Subset.nProf = Input.counts-Subset.pcounts;
            }
            MPI_Send(&Subset.nProf, 1, MPI_INT, pid, TAG_NPROF, \
                MPI_COMM_WORLD);
            MPI_Send(Subset.coord, 2, MPI_INT, pid, TAG_COORD, \
                MPI_COMM_WORLD);

            if(Input.cache_prof){      
              ptr = Input.profbuff+Stk.nw*4*Subset.pcounts;
              PixelMV(&Input, &Subset);
            }else{
              rProfile(&Input, &Stk, &Subset);
            }

            MPI_Send(ptr, Stk.nw*4*Subset.nProf, MPI_DOUBLE, pid, \
                TAG_DATA, MPI_COMM_WORLD);
            pid++;
            actived++;
          }
        }

        if(actived<Mpi.size){
          for(int ii=actived; ii<Mpi.size; ii++){
            MPI_Send(&terminal, 1, MPI_INT, ii, TAG_NPROF, \
                MPI_COMM_WORLD);
          }
        }

        while(rcounts<Input.counts){
          MPI_Recv(&SubsetRV.nProf, 1, MPI_INT, MPI_ANY_SOURCE,
              TAG_R_NPROF, MPI_COMM_WORLD, &status);

          src = status.MPI_SOURCE;
          MPI_Recv(SubsetRV.coord, 2, MPI_INT, src, TAG_R_COORD, \
              MPI_COMM_WORLD, &status);
          if(Input.cache_prof){
            ptrRV = Input.resbuff+10*((SubsetRV.coord[1] \
                -Input.sol_box[1][0])*Input.cache_header.nx \
                +(SubsetRV.coord[0]-Input.sol_box[0][0]));
          }else{
            ptrRV = Input.resbuff;
          }
       
          MPI_Recv(ptrRV, 10*SubsetRV.nProf, MPI_DOUBLE, src, \
              TAG_R_DATA, MPI_COMM_WORLD, &status);

          if(Input.output_fit){
            MPI_Recv(Input.fitbuff, Stk.nw*4*SubsetRV.nProf, MPI_DOUBLE, \
              src, TAG_R_DATA, MPI_COMM_WORLD, &status);  
          }

          if(!Input.cache_prof){
            WRITE_RESULT(&Input, &SubsetRV, &Stk);
          }

          rcounts += SubsetRV.nProf;
          valid = 0;
          while(Subset.pcounts < Input.counts){

            if(!Input.cache_prof && Input.cache[icache]){
              icache++;
              PixelMV(&Input, &Subset);
              rcounts += Subset.nProf;
            }else{
              valid = 1;
              break;
            }
          }

          if(valid){
            icache++;
            MPI_Send(&Subset.nProf, 1, MPI_INT, src, TAG_NPROF, \
                MPI_COMM_WORLD);
            MPI_Send(Subset.coord, 2, MPI_INT, src, TAG_COORD, \
                MPI_COMM_WORLD);

            if(Input.cache_prof){      
              ptr = Input.profbuff+Stk.nw*4*Subset.pcounts;
              PixelMV(&Input, &Subset);
              
            }else{
              rProfile(&Input, &Stk, &Subset);
            }
            MPI_Send(ptr, Stk.nw*4*Subset.nProf, MPI_DOUBLE, \
                src, TAG_DATA, MPI_COMM_WORLD);
      
          }else{
            MPI_Send(&terminal, 1, MPI_INT, src, TAG_NPROF, \
                MPI_COMM_WORLD);
          }
        }
        if(Input.cache_prof){
          Subset.coord[0] = Input.sol_box[0][0];
          Subset.coord[1] = Input.sol_box[1][0];
          Subset.nProf = Input.counts;
          WRITE_RESULT(&Input, &Subset, &Stk);
        }
  
      }else{
        while(1){

          MPI_Recv(&Subset.nProf, 1, MPI_INT, 0, TAG_NPROF, \
              MPI_COMM_WORLD, &status);
          if(Subset.nProf<0) break;

          MPI_Recv(Subset.coord, 2, MPI_INT, 0, TAG_COORD, \
              MPI_COMM_WORLD, &status);
          MPI_Recv(Input.profbuff, Stk.nw*4*Subset.nProf, MPI_DOUBLE, \
              0, TAG_DATA, MPI_COMM_WORLD, &status);  
          INVERSION_MULTI(&Input, &Stk, &Para, &LM, &Subset);
          MPI_Send(&Subset.nProf, 1, MPI_INT, 0, TAG_R_NPROF, \
              MPI_COMM_WORLD);
          MPI_Send(Subset.coord, 2, MPI_INT, 0, TAG_R_COORD, \
              MPI_COMM_WORLD);
          MPI_Send(Input.resbuff, 10*Subset.nProf, MPI_DOUBLE, 0, \
              TAG_R_DATA, MPI_COMM_WORLD);
          if(Input.output_fit){
            MPI_Send(Input.fitbuff, Stk.nw*4*Subset.nProf, MPI_DOUBLE, \
              0, TAG_R_DATA, MPI_COMM_WORLD); 
          }
        }
      }
#endif 
    }

    FREERAM(&Input, &Stk, &LM, &Para, &Mpi);
    Timer(&Mpi);

#ifdef PARALLEL_MPI
    MPI_Finalize();
#endif

    return 0;
}

/*--------------------------------------------------------------------------------*/
