
#include "MPI_INIT.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        06 Mar. 2026  (Hao Li)
          --- Updates:  
              Suport MPI and OpenMP. 
              removed the subroutines Randomseeds and Timer_Mpi, which 
              are redesigned.

        27 Nov. 2024
          --- Updates:  a new subroutine Randomseeds to generate random 
                        seeds for mult-processor (Hao Li)
                        a new subroutine Timer_Mpi to compute and 
                        print the running time (Hao Li)

        27 Apr. 2024
          --- Initial commit (Hao Li)
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

int cpu_check(int *cpu_busy, int ncpu){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        return the id of a free cpu if there is.
      Record of revisions:
        25 Apr. 2024
      Input parameters:
        cpu_busy, an array with the cpu status.
        ncpu, the length of the array.
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/
    
    int icpu;

    for(icpu=1;icpu<ncpu;icpu++){
      if(!cpu_busy[icpu]) return icpu;
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

void MPI_SETUP(STRUCT_MPI *Mpi){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        set up the rank and mpi size.
      Record of revisions:
        02 Mar. 2026.
      Input parameters:
        Mpi, a structure storing mpi information.
      Output parameters:
        Mpi, a structure storing mpi information.        
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

#ifdef USE_MPI
    Mpi->comm = MPI_COMM_WORLD;
    MPI_Comm_rank(Mpi->comm, &Mpi->rank);
    MPI_Comm_size(Mpi->comm, &Mpi->size);
#else
    Mpi->rank = 0;
    Mpi->size = 1;
#endif

    return;
}

/*--------------------------------------------------------------------------------*/

void OPENMP_SETUP(STRUCT_MPI *Mpi){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        set up the thread id and numbers.
      Record of revisions:
        02 Mar. 2026.
      Input parameters:
        Mpi, a structure storing mpi info.
      Output parameters:
        Mpi, a structure storing mpi info.        
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

#ifdef USE_OPENMP
    omp_set_num_threads(Mpi.nthreads);
    Mpi->thread_id = omp_get_thread_num(); 
#else 
    Mpi->thread_id = 0; 
#endif

    return;
}

/*--------------------------------------------------------------------------------*/
