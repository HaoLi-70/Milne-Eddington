
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <mpi.h>
#include <fitsio.h>
#include "ALLOCATION.h"
#include "MILNE_EDDINGTON.h"
#include "LM_FIT.h"
#include "READ_INPUT.h"
#include "MPI_CTRL.h"
#include "TIME_PRINT.h"
#include "IO.h"
#include "RANDOM_NUMBER.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        27 Nov. 2024
          --- Updates:  a new subroutine to generate random seeds for 
                        mult-processor (Hao Li)
                        fast mode of the inversion (Hao Li)

        28 Jun. 2024
          --- Initial commit (Hao Li)
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

#define Input_Path "./input"

/*--------------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Stokes inversion (Milne-Eddington).
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    MPI_Init(&argc, &argv);
    STRUCT_MPI *Mpi = (STRUCT_MPI *)malloc(sizeof(STRUCT_MPI));
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &(Mpi->rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(Mpi->nprocs));

    Time_Print_Mpi(Mpi);

    STRUCT_INPUT *Input = (STRUCT_INPUT *)malloc(sizeof(STRUCT_INPUT));

    char *Filename = Input_Path;
    if(argc>=2){
      Filename = argv[1];
    }

    if(!FILE_EXIST(Filename)){
      free(Input);
      if (Mpi->rank == 0){
        Error(enum_error, "main", "input doesn't exist! \n", NULL);
      }else{
        ABORTED();
      }
    }

    Mrank = Mpi->rank;
    Mnprocs = Mpi->nprocs;
    int *cpu_busy = (int *)VECTOR(1, Mpi->nprocs-1, enum_int, true);

    STRUCT_PAR *Par = (STRUCT_PAR *)malloc(sizeof(STRUCT_PAR));
    STRUCT_STK *Stk = (STRUCT_STK *)malloc(sizeof(STRUCT_STK));
    STRUCT_LM *LM = (STRUCT_LM *)malloc(sizeof(STRUCT_LM));

    // read the input configuration.
    RDINPUT(Filename, Input, Stk, LM, Par, Mpi);

    // read the wavelength points. 
    rWavelength(Input, Stk, Mpi);

    long fpixel[4], naxes[4];
    int coord[2], coordr[2], ipid, scounts = 0, rcounts = 0, ecounts = 0;
    int i, ix, iy, receiving, stat = 0, naxis, iw, ipar;
    double Imean, ***res;
    float **profall;
    char filename[Max_Line_Length+1];

    Randomseeds(Mpi);

    if(Input->fastmode){

      Input->nxcache = 1+Input->sol_box[0][1]-Input->sol_box[0][0];
      Input->nycache = 1+Input->sol_box[1][1]-Input->sol_box[1][0];
      Input->counts = Input->nxcache*Input->nycache;

      if(Mpi->rank==0){
        Stk->profall = (float ***)TENSOR_FLT(0, Input->ny-1, \
              0, Input->nx-1, 0, 4*Stk->nl-1, false);

        rprofileall(Input, Stk);

        res = (double ***)TENSOR_DBL(0, Input->nycache-1, \
                0, Input->nxcache-1, 0, 9, true);
      }else{

        profall = (float **)MATRIX(0, 3, 0, Stk->nl-1, enum_flt, false);

      }

      if(Mpi->nprocs>1){

        if(Mpi->rank==0){

          coord[0] = Input->sol_box[0][0];
          coord[1] = Input->sol_box[1][0];
          while(true){

            /* if aborting */
            //if(true){


            //}

            /* code */
            // if there are pixel left and free CPU and 
            ipid = cpu_check(cpu_busy, Mpi->nprocs);
      
            if(scounts<Input->counts && ipid>0){

              ix = coord[0]-Input->sol_box[0][0];
              iy = coord[1]-Input->sol_box[1][0];

              Imean = 0.0;
              for(iw=0;iw<Stk->nl;iw++){
                Imean += Stk->profall[coord[1]][coord[0]][iw];
              }
              Imean /= Stk->nl;

              // if on disk
              if(Imean>Input->Icriteria){
                // send the coordinates
                MPI_Send(coord, 2, MPI_INT, ipid, 1000+ipid, MPI_COMM_WORLD);

                // send the data
                MPI_Send(Stk->profall[coord[1]][coord[0]], Stk->nl*4, \
                    MPI_FLOAT, ipid, 1000+ipid, MPI_COMM_WORLD);

                // That CPU is now busy
                cpu_busy[ipid] = 1;

              // if off limb
              }else{

                rcounts++;

              }
        
              // move to the next pixel
              if(coord[0]>=Input->sol_box[0][1]){
                coord[0] = Input->sol_box[0][0];
                coord[1]++;
              }else{
                coord[0]++;
              }

              scounts++;

            }

            for (ipid = 1; ipid < Mpi->nprocs; ipid++){

              //Test for a message from the slave
              MPI_Iprobe(ipid, 2000+ipid, MPI_COMM_WORLD, &receiving, \
                  &status);
              if(!receiving) continue;

              // recieve the coordinates
              MPI_Recv(coordr, 2, MPI_INT, ipid, 2000+ipid, MPI_COMM_WORLD, \
                  &status);

              // coordiante of the result array
              coordr[0] = coordr[0]-Input->sol_box[0][0];
              coordr[1] = coordr[1]-Input->sol_box[1][0];

              // recieve the result
              MPI_Recv(res[coordr[1]][coordr[0]], 9, MPI_DOUBLE, ipid, \
                  2000+ipid, MPI_COMM_WORLD, &status);

              MPI_Recv(res[coordr[1]][coordr[0]]+9, 1, MPI_DOUBLE, ipid, \
                  2000+ipid, MPI_COMM_WORLD, &status);

              rcounts++;
              // the cpu is free
              cpu_busy[ipid] = 0;
            }

            if(rcounts>=Input->counts) break;

          }
        
          coord[0] = -1;
          for (ipid = 1; ipid < Mpi->nprocs; ipid++){
            // notice the slave
            MPI_Send(coord, 2, MPI_INT, ipid, 1000+ipid, MPI_COMM_WORLD);
          }

          sprintf(filename, "!%s",Input->Result_Path);
          fits_create_file(&(Input->fptr_res), filename, &stat);

          naxis = 3;
          naxes[0] = 10;
          naxes[1] = Input->nxcache;
          naxes[2] = Input->nycache;

          fits_create_img(Input->fptr_res, DOUBLE_IMG, naxis, naxes, &stat);

          //write results
          fits_open_file(&(Input->fptr_res), Input->Result_Path,\
              READWRITE, &stat);

          fpixel[0] = 1;
          fpixel[1] = 1;
          fpixel[2] = 1;
              
          fits_write_pix(Input->fptr_res, TDOUBLE, fpixel, \
              10*Input->nxcache*Input->nycache, res[0][0], &stat);;
              
          fits_close_file(Input->fptr_res, &stat);


        }else{

          while(true){
            // recieve the integers
            MPI_Recv(coord, 2, MPI_INT, 0, 1000+Mpi->rank, MPI_COMM_WORLD, \
                &status);
            // break if no data
            if(coord[0]<0) break;

            // recieve the data
            MPI_Recv(profall[0], Stk->nl*4, MPI_FLOAT, 0, 1000+Mpi->rank,\
                MPI_COMM_WORLD, &status); 
            
            for(iw=0;iw<Stk->nl;iw++){
              Stk->prof[0][iw] = profall[0][iw];
              Stk->prof[1][iw] = profall[1][iw];
              Stk->prof[2][iw] = profall[2][iw];
              Stk->prof[3][iw] = profall[3][iw];
            }

            INVERSION(Input, Stk, Par, LM, Mpi); 

            MPI_Send(coord, 2, MPI_INT, 0, 2000+Mpi->rank, MPI_COMM_WORLD);

            MPI_Send(Par->Par_Best+1, 9, MPI_DOUBLE, 0, 2000+Mpi->rank, \
                MPI_COMM_WORLD);

            MPI_Send(&Par->Chisq_Best, 1, MPI_DOUBLE, 0, 2000+Mpi->rank, \
                MPI_COMM_WORLD);

          }
        }

      }else{

        for(iy=0;iy<Input->nycache;iy++){
          if(iy%50 == 0){
            sprintf(MeSS, " -- pixel y = %d -- \n", iy);
            VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
          }
      
          for(ix=0;ix<Input->nxcache;ix++){

            coord[0] = ix+Input->sol_box[0][0];
            coord[1] = iy+Input->sol_box[1][0];

            Imean = 0.0;
            for(iw=0;iw<Stk->nl;iw++){
              Imean += Stk->profall[coord[1]][coord[0]][iw];
            }
            Imean /= Stk->nl;

            if(Imean!=Imean) continue;

            // if on disk
            if(Imean<Input->Icriteria) continue;

            for(iw=0;iw<Stk->nl;iw++){
              Stk->prof[0][iw] = (double) \
                  (Stk->profall[coord[1]][coord[0]][iw]);
              Stk->prof[1][iw] = (double) \
                  (Stk->profall[coord[1]][coord[0]][iw+Stk->nl]);
              Stk->prof[2][iw] = (double) \
                  (Stk->profall[coord[1]][coord[0]][iw+Stk->nl*2]);
              Stk->prof[3][iw] = (double) \
                  (Stk->profall[coord[1]][coord[0]][iw+Stk->nl*3]);
            }

            INVERSION(Input, Stk, Par, LM, Mpi); 

            for(ipar=0;ipar<9;ipar++){
              res[iy][ix][ipar] = Par->Par_Best[ipar+1];
            }
            res[iy][ix][9] = Par->Chisq_Best;
   
          }
        }

        sprintf(filename, "!%s",Input->Result_Path);
        fits_create_file(&(Input->fptr_res), filename, &stat);

        naxis = 3;
        naxes[0] = 10;
        naxes[1] = Input->nxcache;
        naxes[2] = Input->nycache;

        fits_create_img(Input->fptr_res, DOUBLE_IMG, naxis, naxes, &stat);

        //write results
        fits_open_file(&(Input->fptr_res), Input->Result_Path,\
            READWRITE, &stat);

        fpixel[0] = 1;
        fpixel[1] = 1;
        fpixel[2] = 1;
              
        fits_write_pix(Input->fptr_res, TDOUBLE, fpixel, \
            10*Input->nxcache*Input->nycache, res[0][0], &stat);
              
        fits_close_file(Input->fptr_res, &stat);

      }

    }else{

      int **cache = cache_init(Input, Mpi, Stk, &stat);
      
      MPI_Bcast(&stat, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stat!=0){
        if(Mpi->rank==0){
          FREE_MATRIX(cache, 0, 0, enum_int);  
        }

        Free_Ram(LM, Stk, Input, Par, Mpi);

        if(Mpi->rank==0){
          Error(enum_error, "main", "error return from cache_init \n", \
              Input->Verbose_Path);
        }else{
          ABORTED();
        }
      }

      if(Mpi->rank==0){
        sprintf(MeSS, "\n -- inversion starting with %d cores -- \n", \
            Mpi->nprocs);
        VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
      }

      if(Mpi->nprocs>1){

        if(Mpi->rank==0){

          coord[0] = Input->sol_box[0][0];
          coord[1] = Input->sol_box[1][0];
          while(true){

            /* if aborting */
            //if(true){


            //}

            /* code */
            // if there are pixel left and free CPU and 
            ipid = cpu_check(cpu_busy, Mpi->nprocs);
      
            if(scounts<Input->counts && ipid>0){
              ix = coord[0]-Input->sol_box[0][0];
              iy = coord[1]-Input->sol_box[1][0];

              if(cache[ix][iy] == 0){

                // read the data 
                rprofile(Input, coord, Stk);

                Imean = 0.0;
                for(i=0;i<Stk->nl;i++) Imean += Stk->prof[0][i];
                Imean /= Stk->nl;

                // if on disk
                if(Imean>Input->Icriteria){
                  // send the coordinates
                  MPI_Send(coord, 2, MPI_INT, ipid, 1000+ipid, MPI_COMM_WORLD);

                  // send the data
                  MPI_Send(Stk->prof[0], Stk->nl*4, MPI_DOUBLE, ipid, \
                      1000+ipid, MPI_COMM_WORLD);

                  // That CPU is now busy
                  cpu_busy[ipid] = 1;

                // if off limb
                }else{

                  //write cache
                  cache_write(Input, coordr);

                  rcounts++;
                }

              }else{
                rcounts++;
              }
                          
              // move to the next pixel
              if(coord[0]>=Input->sol_box[0][1]){
                coord[0] = Input->sol_box[0][0];
                coord[1]++;
              }else{
                coord[0]++;
              }

              scounts++;

            }

            for (ipid = 1; ipid < Mpi->nprocs; ipid++){

              //Test for a message from the slave
              MPI_Iprobe(ipid, 2000+ipid, MPI_COMM_WORLD, &receiving, \
                  &status);
              if(!receiving) continue;

              // recieve the coordinates
              MPI_Recv(coordr, 2, MPI_INT, ipid, 2000+ipid, MPI_COMM_WORLD, \
                  &status);

              coordr[0] = coordr[0]-Input->sol_box[0][0];
              coordr[1] = coordr[1]-Input->sol_box[1][0];

              // recieve the result
              MPI_Recv(Par->Par_Best+1, 9, MPI_DOUBLE, ipid, 2000+ipid, \
                  MPI_COMM_WORLD, &status);

              MPI_Recv(&Par->Chisq_Best, 1, MPI_DOUBLE, ipid, 2000+ipid, \
                  MPI_COMM_WORLD, &status);

              //write results
              fits_open_file(&(Input->fptr_res), Input->Result_Path,\
                  READWRITE, &stat);

              fpixel[1] = coordr[0]+1;
              fpixel[2] = coordr[1]+1;
              fpixel[0] = 1;
              fits_write_pix(Input->fptr_res, TDOUBLE, fpixel, 9, \
                  Par->Par_Best+1, &stat);

              fpixel[0] = 10;
              fits_write_pix(Input->fptr_res, TDOUBLE, fpixel, 1, \
                  &Par->Chisq_Best, &stat);
              
              fits_close_file(Input->fptr_res, &stat);

              if(Input->output_fit){

                // recieve the fit
                MPI_Recv(Stk->syn[0], Stk->nl*4, MPI_DOUBLE, ipid, \
                    2000+ipid, MPI_COMM_WORLD, &status);

                // write the fit
                fits_open_file(&(Input->fptr_fit), Input->Fit_Path, \
                    READWRITE, &stat);

                // set the position of the first pixel.
                fpixel[0] = 1;
                fpixel[1] = 1;
                fpixel[2] = coordr[0]+1;
                fpixel[3] = coordr[1]+1;
                fits_write_pix(Input->fptr_fit, TDOUBLE, fpixel, Stk->nl*4, \
                    Stk->syn[0], &stat);

                fits_close_file(Input->fptr_fit, &stat);

              }   

              //write cache
              cache_write(Input, coordr);


              rcounts++;
              // the cpu is free
              cpu_busy[ipid] = 0;
            }

            if(rcounts>=Input->counts) break;

          }
        
          coord[0] = -1;
          for (ipid = 1; ipid < Mpi->nprocs; ipid++){
              // notice the slave
              MPI_Send(coord, 2, MPI_INT, ipid, 1000+ipid, MPI_COMM_WORLD);
          }

          FREE_MATRIX(cache, 0, 0, enum_int);  

        }else{

          while(true){
            // recieve the integers
            MPI_Recv(coord, 2, MPI_INT, 0, 1000+Mpi->rank, MPI_COMM_WORLD, \
                &status);
            // break if no data
            if(coord[0]<0) break;

            // recieve the data
            MPI_Recv(Stk->prof[0], Stk->nl*4, MPI_DOUBLE, 0, 1000+Mpi->rank,\
                MPI_COMM_WORLD, &status);  

            INVERSION(Input, Stk, Par, LM, Mpi); 
            MPI_Send(coord, 2, MPI_INT, 0, 2000+Mpi->rank, MPI_COMM_WORLD);

            MPI_Send(Par->Par_Best+1, 9, MPI_DOUBLE, 0, 2000+Mpi->rank, \
                MPI_COMM_WORLD);

            MPI_Send(&Par->Chisq_Best, 1, MPI_DOUBLE, 0, 2000+Mpi->rank, \
                MPI_COMM_WORLD);
            if(Input->output_fit){
              MPI_Send(Stk->syn[0], Stk->nl*4, MPI_DOUBLE, 0, \
                  2000+Mpi->rank, MPI_COMM_WORLD);

            }
          }
        }

      }else{

        for(coordr[1]=0;coordr[1]<Input->nycache;coordr[1]++){
          if(coordr[1]%50 == 0){
            sprintf(MeSS, " -- pixel y = %d -- \n", coordr[1]);
            VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
          }
      
          for(coordr[0]=0;coordr[0]<Input->nxcache;coordr[0]++){
            coord[0] = coordr[0]+Input->sol_box[0][0];
            coord[1] = coordr[1]+Input->sol_box[1][0];
  /*
            if(coordr[0]%50 == 0){

              sprintf(MeSS, " -- pixel x = %d -- \n", coordr[0]);
              VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>3);
            }
  */
            if(cache[coordr[0]][coordr[1]] == 0){
              // read the data 
              coord[0] = coordr[0]+Input->sol_box[0][0];
              coord[1] = coordr[1]+Input->sol_box[1][0];
              rprofile(Input, coord, Stk);

              Imean = 0.0;
              for(i=0;i<Stk->nl;i++) Imean += Stk->prof[0][i];
              Imean /= Stk->nl;

              // if on disk
              if(Imean<Input->Icriteria){

                //write cache
                cache_write(Input, coordr);

                continue;
              }

              INVERSION(Input, Stk, Par, LM, Mpi); 

              fits_open_file(&(Input->fptr_res), Input->Result_Path, \
                  READWRITE, &stat);

              fpixel[1] = coordr[0]+1;
              fpixel[2] = coordr[1]+1;

              fpixel[0] = 1;
              fits_write_pix(Input->fptr_res, TDOUBLE, fpixel, 9, \
                  Par->Par_Best+1, &stat);

              fpixel[0] = 10;

              fits_write_pix(Input->fptr_res, TDOUBLE, fpixel, 1, \
                  &Par->Chisq_Best, &stat);       

              fits_close_file(Input->fptr_res, &stat);

              
              cache_write(Input, coordr);

              if(Input->output_fit){

                // write the fit
                fits_open_file(&(Input->fptr_fit), Input->Fit_Path, \
                    READWRITE, &stat);
                fpixel[0] = 1;
                fpixel[1] = 1;
                fpixel[2] = coordr[0]+1;
                fpixel[3] = coordr[1]+1;
                fits_write_pix(Input->fptr_fit, TDOUBLE, fpixel, Stk->nl*4, \
                    Stk->syn[0], &stat);

                fits_close_file(Input->fptr_fit, &stat);

              }   
          
            }
          }
        }

        FREE_MATRIX(cache, 0, 0, enum_int);  

      }
    }

    Time_Print_Mpi(Mpi);

    if (Mpi->rank == 0){ 
      sprintf(MeSS, "\n -- inversion finished -- \n");
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }

    //Free_Ram(LM, Stk, Input, Par, Mpi);

    MPI_Finalize();

    return 0;
    
}

/*--------------------------------------------------------------------------------*/
