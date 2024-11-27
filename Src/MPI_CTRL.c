
#include "MPI_CTRL.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        27 Nov. 2024
          --- Updates:  a new subroutine Randomseeds to generate random 
                        seeds for mult-processor (Hao Li)
                        a new subroutine Time_Print_Mpi to compute and 
                        print the running time (Hao Li)

        27 Apr. 2024
          --- Initial commit (Hao Li)
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

int Mrank;
int Mnprocs;

extern int cpu_check(int *cpu_busy, int ncpu){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        return the id of a free cpu if there is.
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        cpu_busy, an array with the cpu status.
        ncpu, the length of the array.
      Return:
        .
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
    
    // true = 1, flase = 0;
    int icpu;

    for(icpu=1;icpu<ncpu;icpu++){
      if(!cpu_busy[icpu]) return icpu;
    }

    return 0;

}

/*--------------------------------------------------------------------------------*/

extern void CONTROL(void){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        controls if any CPU has crashed and stops if needed..
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        .
      Return:
        .
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    int a1=0, a2=0;
    
    MPI_Allreduce(&a1, &a2, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    
    if (!a2) return;
    
    MPI_Finalize();
    
    return;
}

/*--------------------------------------------------------------------------------*/

extern void ABORTED(void){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        terminate the job.
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        .
      Return:
        .
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
  
  MPI_Finalize();

  exit(0);

}

/*--------------------------------------------------------------------------------*/

extern int Randomseeds(STRUCT_MPI *Mpi){

    /*######################################################################
      Purpose:
        initialize the mpi structure.
      Record of revisions:
        30 Otc. 2024
      Input parameters:
        Mpi, a structure save the information for mpi.
      Output parameters:
        Mpi, the structure save the information for mpi.
    ######################################################################*/
    
    int ipid;

    Mpi->idum = (long *)malloc(sizeof(long));
    long *tmp = (long *)malloc(sizeof(long)*Mpi->nprocs);

    if (Mpi->rank == 0){
      srand((unsigned int)time(NULL));
      srand(rand());  
      for(ipid=0; ipid<Mpi->nprocs; ipid++){  
        tmp[ipid] = -rand();
      }
    }
    MPI_Bcast(tmp, Mpi->nprocs, MPI_LONG, 0, MPI_COMM_WORLD);

    *(Mpi->idum) = tmp[Mpi->rank];
    free(tmp);

    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int Time_Print_Mpi(STRUCT_MPI *Mpi){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        compute and print the running time.
      Record of revisions:
        27 Nov. 2024 (Hao Li)
      Input parameters:
        Mpi, structure with MPI information.
      Return:
        return the current conts.
    ######################################################################*/
    
/*--------------------------------------------------------------------------------*/

    static int conts = 0;
    static clock_t time_begin = 0;
    
    if(conts == 0){
      time_begin = clock();
      if(Mpi->rank == 0){
        fprintf(stderr, "\n Time Calculating Is Initialized \n");
      }
      conts++;
        
    }else{
      long time_tmp, timec = clock()-time_begin;

      MPI_Reduce(&timec, &time_tmp, 1, MPI_LONG, MPI_MAX, 0, \
          MPI_COMM_WORLD);

      long hours = time_tmp/CLOCKS_PER_SEC/3600;      
      long minutes = (time_tmp/CLOCKS_PER_SEC-hours*3600)/60;
      double seconds = (double)(time_tmp)/CLOCKS_PER_SEC \
          -hours*3600-minutes*60;
        
      if(Mpi->rank == 0){
        fprintf(stderr, "\n Time Print Point %d \n",conts);
          
        if(hours > 0){
          fprintf(stderr," Running time= %lu h %lu min %.2lf sec\n", \
              hours,minutes,seconds);
        }else if(minutes > 0){
          fprintf(stderr," Running time= %lu min %.2lf sec\n", \
              minutes, seconds);
        }else{
          fprintf(stderr," Running time= %.2lf sec\n", seconds);
        }
      }
      conts++;
    }
    return conts-1;
}

/*--------------------------------------------------------------------------------*/