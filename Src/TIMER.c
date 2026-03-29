
#include "TIMER.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        15 Mar. 2026  (Hao Li)
          --- Updates:  
              Suport MPI and OpenMP. 

        8 Sept. 2021
          --- Initial commit (Hao Li)
     
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

int Timer(STRUCT_MPI *Mpi){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        compute and display the execution time.
      Record of revisions:
        15 Mar. 2026. 
      Input parameters:
         Mpi, a structure storing Mpi info.
      Return:
        return the current conts.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    static int count = 0;
    static double t_begin = 0.0;
    double t_global = 0; 

    if(count == 0){

#if defined(USE_MPI)
      t_begin = MPI_Wtime();
#elif defined(UST_OPENMP)
      t_global = omp_get_wtime() - t_begin;
#else
      t_begin = (double)clock()/CLOCKS_PER_SEC;
#endif

      if(Mpi->rank == 0) fprintf(stderr,"\n Timer initialized\n");
      count++;

    }else{

#if defined(USE_MPI)
      double t_local = MPI_Wtime() - t_begin;
      MPI_Reduce(&t_local, &t_global, 1, MPI_DOUBLE, MPI_MAX,
          0,MPI_COMM_WORLD);
#elif defined(UST_OPENMP)
      t_global = omp_get_wtime() - t_begin;
#else
      t_global = (double)clock()/CLOCKS_PER_SEC-t_begin;
#endif

      if(Mpi->rank == 0){

        long hours = (long)(t_global/3600);
        long minutes = (long)((t_global-hours*3600)/60);
        double sec = t_global-hours*3600-minutes*60;

        fprintf(stderr,"\n Timer display %d\n",count);

        if(hours > 0){
          fprintf(stderr, " Execution time = %ld h %ld min %.2f sec\n",
              hours, minutes, sec);
        }else if(minutes > 0){
          fprintf(stderr, " Execution time = %ld min %.2f sec\n",
              minutes, sec);
        }else{
          fprintf(stderr, " Execution time = %.2f sec\n", sec);
        }
      }

      count++;
    }

    return count - 1;
}

/*--------------------------------------------------------------------------------*/

