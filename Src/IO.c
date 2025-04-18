
#include "IO.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        17 Apr. 2025
          --- Updates: locate the line center (Hao Li)

        11 Apr. 2025
          --- bugfix:  does not save the best fit (Hao Li)

        27 Nov. 2024
          --- Updates: a new subroutine rprofileall to read all the 
              profiles for the fast mode (Hao Li)

        28 Jun. 2024
          --- Initial commit (Hao Li)

     ######################################################################*/

/*--------------------------------------------------------------------------------*/

extern int rWavelength(STRUCT_INPUT *Input, STRUCT_STK *Stk, \
        STRUCT_MPI *Mpi){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        read the wavelength.
      Record of revisions:
        17 Apr. 2025 (Hao Li)
      Input parameters:
        Input, the input configuration.
        Stk, structure with Stokes profiles.
        Mpi, structure with mpi configuration.
      Output parameters:
        Stk, structure with Stokes profiles.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    const char *rname = "rWavelength";

    // CFITSIO status value MUST be initialized to zero!
    int status = 0, i, j;

    // number of the HDUs, the bit of each pixel, the number of axes
    int bitpix, naxis;
    
    // the length of each axis, the first pixel
    long naxes[4], fpixel = 1;

    float *tmp;

    double Imean;
    
    if(Mpi->rank==0){

      fits_open_file(&(Input->fptr_data), Input->Data_Path, \
          READONLY, &status);

      if(status!=0){
        sprintf(MeSS, "error in opening the data file.\n");
        goto ERR1;        
      }

      fits_open_file(&(Input->fptr_wav), Input->Wav_Path, \
          READONLY, &status);

      if(status!=0){
        sprintf(MeSS, "error in opening the wavelength file.\n");
        goto ERR1;        
      }
      
      // get the dimension of the data hdu
      fits_get_img_dim(Input->fptr_data, &(Input->naxis), &status);

      if(status!=0){
        sprintf(MeSS, "error in getting the dimension of the data hdu.\n");
        goto ERR1;        
      }

      if(Input->naxis<3||Input->naxis>4){ 
        sprintf(MeSS, "error in the dimension of the data hdu.\n");
        goto ERR1;      
      }

      // get the size of each dimension 
      fits_get_img_size(Input->fptr_data, Input->naxis, naxes, &status);

      if(status!=0){
        sprintf(MeSS, "error in getting the size of the data hdu.\n");
        goto ERR1;        
      }

      if(naxes[1]!=4){
        sprintf(MeSS, "error in the size of the data hdu.\n");
        goto ERR1;              
      }

      Stk->nl = naxes[0];
      Input->nx = naxes[2];
      if(Input->naxis==3){
        Input->ny = 1;        
      }else{
        Input->ny = naxes[3];
      }

      // check the solution box
      if(Input->sol_box[0][0]<0||Input->sol_box[0][0]>=Input->nx) \
          Input->sol_box[0][0] = 0;
      if(Input->sol_box[0][1]<0||Input->sol_box[0][1]>=Input->nx) \
          Input->sol_box[0][1] = Input->nx-1;
      if(Input->sol_box[1][0]<0||Input->sol_box[1][0]>=Input->ny) \
          Input->sol_box[1][0] = 0;
      if(Input->sol_box[1][1]<0||Input->sol_box[1][1]>=Input->ny) \
          Input->sol_box[1][1] = Input->ny-1;

      // get the type of the hdu
      fits_get_img_type(Input->fptr_data, &bitpix, &status);

      if(status!=0){
        sprintf(MeSS, "error in getting the data type of data hdu.\n");
        goto ERR1;        
      }

      Stk->prof = (double **)MATRIX(0, 3, 0, Stk->nl-1, enum_dbl, false);
      
      //Input->Icriteria = 0.;
      if(bitpix==16){ 
        Input->type_data = enum_int;  
        Input->data_int = (int *)malloc(4*Stk->nl*sizeof(int));
        /*
        if(Input->nx>=Input->ny){
          if(Input->ny>1){
            Input->fpixel[3] = Input->ny/2;
          }else{
            Input->fpixel[3] = 1;
          }
          for(i=1;i<=Input->nx;i++){
            Input->fpixel[2] = i;
            fits_read_pix(Input->fptr_data, TSHORT, Input->fpixel, \
                Stk->nl, NULL, Input->data_int, NULL, &status); 
            Imean = 0.0;
            for(j=0;j<Stk->nl;j++){
               Imean += Input->data_int[j];
            }
            Imean /= Stk->nl;
            if(Input->Icriteria<Imean) Input->Icriteria = Imean;
          }
          
        }else{
          if(Input->nx>1){
            Input->fpixel[2] = Input->nx/2;
          }else{
            Input->fpixel[2] = 1;
          }
          for(i=1;i<=Input->ny;i++){
            Input->fpixel[3] = i;
            fits_read_pix(Input->fptr_data, TSHORT, Input->fpixel, \
                Stk->nl, NULL, Input->data_int, NULL, &status); 
            for(j=0;j<Stk->nl;j++){
              if(Input->Icriteria<Input->data_int[j]){ 
                Input->Icriteria = (double)Input->data_int[j];
              }
            }
          }
        }
        */
      }else if(bitpix==-32){
        Input->type_data = enum_flt;
        Input->data_flt = (float *)malloc(4*Stk->nl*sizeof(float));
        /*
        if(Input->nx>=Input->ny){
          if(Input->ny>1){
            Input->fpixel[3] = Input->ny/2;
          }else{
            Input->fpixel[3] = 1;
          }
          for(i=1;i<=Input->nx;i++){
            Input->fpixel[2] = i;
            fits_read_pix(Input->fptr_data, TFLOAT, Input->fpixel, \
                Stk->nl, NULL, Input->data_flt, NULL, &status);

            Imean = 0.0;
            for(j=0;j<Stk->nl;j++) Imean += Input->data_flt[j];
            Imean /= Stk->nl;
            if(Input->Icriteria<Imean) Input->Icriteria = Imean;
            
          }

        }else{
          if(Input->nx>1){
            Input->fpixel[2] = Input->nx/2;
          }else{
            Input->fpixel[2] = 1;
          }
          for(i=1;i<=Input->ny;i++){
            Input->fpixel[3] = i;
            fits_read_pix(Input->fptr_data, TFLOAT, Input->fpixel, \
                Stk->nl, NULL, Input->data_flt, NULL, &status); 
            for(j=0;j<Stk->nl;j++){
              if(Input->Icriteria<Input->data_flt[j]){ 
                Input->Icriteria = (double)Input->data_flt[j];
              }
            }
          }
        }
        */
      }else if(bitpix==-64){
      
        Input->type_data = enum_dbl;
        /*
        if(Input->nx>=Input->ny){
          if(Input->ny>1){
            Input->fpixel[3] = Input->ny/2;
          }else{
            Input->fpixel[3] = 1;
          }

          for(i=1;i<=Input->nx;i++){
            Input->fpixel[2] = i;

            fits_read_pix(Input->fptr_data, TDOUBLE, Input->fpixel, \
                Stk->nl, NULL, Stk->prof[0], NULL, &status); 

            Imean = 0.0;
            for(j=0;j<Stk->nl;j++) Imean += Stk->prof[0][j];
            Imean /= Stk->nl;
            if(Input->Icriteria<Imean) Input->Icriteria = Imean;
            
          }

        }else{
          if(Input->nx>1){
            Input->fpixel[2] = Input->nx/2;
          }else{
            Input->fpixel[2] = 1;
          }
          for(i=1;i<=Input->ny;i++){
            Input->fpixel[3] = i;
            fits_read_pix(Input->fptr_data, TDOUBLE, Input->fpixel, \
                Stk->nl, NULL, Stk->prof[0], NULL, &status);
            for(j=0;j<Stk->nl;j++){
              if(Input->Icriteria<Stk->prof[0][j]){ 
                Input->Icriteria = (double)Stk->prof[0][j];
              }
            }
          }
        }
        */
      }else{
        sprintf(MeSS, "error in the data type of data hdu.\n");
        goto ERR1;  
      }

      //Input->Icriteria *= 0.05;
      

      // get the dimension of the wavlength hdu
      fits_get_img_dim(Input->fptr_wav, &naxis, &status);
      if(status!=0){
        sprintf(MeSS, "error in getting the dimension of wavelength hdu.\n");
        goto ERR1;        
      }

      if(naxis!=1){ 
        sprintf(MeSS, "error in the dimension of wavelength hdu.\n");
        goto ERR1;   
      }

      // get the size of the wavelength 
      fits_get_img_size(Input->fptr_wav, naxis, naxes, &status);
      if(status!=0){
        sprintf(MeSS, "error in getting the size of the wavelength.\n");
        goto ERR1;        
      }

      if(Stk->nl!=naxes[0]){
        sprintf(MeSS, "wavelength is not consistent with the profile.\n");
        goto ERR1;              
      }


      // get the type of the hdu
      fits_get_img_type(Input->fptr_wav, &bitpix, &status);

      if(status!=0){
        sprintf(MeSS, "error in getting the data type of the wavelength hdu.\n");
        goto ERR1;        
      }

      if(bitpix==-32){
        Input->type_wav = enum_flt;
      }else if(bitpix==-64){
        Input->type_wav = enum_dbl;
      }else{
        sprintf(MeSS, "error in the data type of the wavelength hdu.\n");
        goto ERR1;  
      }

      Stk->Lambda = (double *)malloc(Stk->nl*sizeof(double));

      // read wavelength
      if(Input->type_wav == enum_dbl){
        fits_read_pix(Input->fptr_wav, TDOUBLE, &fpixel, Stk->nl, NULL, \
            Stk->Lambda, NULL, &status);        
      }else{
        tmp = (float *)malloc(Stk->nl*sizeof(float));
        fits_read_pix(Input->fptr_wav, TFLOAT, &fpixel, Stk->nl, NULL, \
            tmp, NULL, &status);
        for(i=0;i<Stk->nl;i++) Stk->Lambda[i] = (double)tmp[i];
        free(tmp);
      }

      if(status!=0){
        sprintf(MeSS, "error in reading the wavelength.\n");
        goto ERR1;        
      }

      fits_close_file(Input->fptr_data, &status);

      if(status!=0){
        sprintf(MeSS, "error in closing the data file.\n");
        goto ERR1;        
      }

      fits_close_file(Input->fptr_wav, &status);

      if(status!=0){
        sprintf(MeSS, "error in closing the wavelength file.\n");
        goto ERR1;        
      }

      // broadcast to slaves
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&Stk->nl, 1, MPI_INT, 0, MPI_COMM_WORLD);

      // broadcast the wavelengths
      MPI_Bcast(Stk->Lambda, Stk->nl, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }else{

      // if error in reading wavelength
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(status!=0) ABORTED();

      //  receive the size of the wavelength
      MPI_Bcast(&Stk->nl, 1, MPI_INT, 0, MPI_COMM_WORLD);
      Stk->Lambda = (double *)malloc(Stk->nl*sizeof(double));

      // receive the wavelength
      MPI_Bcast(Stk->Lambda, Stk->nl, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      Stk->prof = (double **)MATRIX(0, 3, 0, Stk->nl-1, enum_dbl, false);
    }

    for(i=0;i<Stk->nl-1;i++){
      if(Input->Lines[0].Lambda0>=Stk->Lambda[i] && 
          Input->Lines[0].Lambda0<=Stk->Lambda[i+1]){
        if(Stk->Lambda[i+1]-Input->Lines[0].Lambda0 <= 
            Input->Lines[0].Lambda0-Stk->Lambda[i]){
          Input->iw = i;
        }else{
          Input->iw = i+1;
        }
      }
    }

    Stk->syn = (double **)MATRIX(0, 3, 0, Stk->nl-1, enum_dbl, false);
    Stk->fit = (double **)MATRIX(0, 3, 0, Stk->nl-1, enum_dbl, false);

    if(Mpi->nprocs==1){
      Stk->Jacobian = (double ***)TENSOR_DBL(0, 3, 1, \
          9, 0, Stk->nl-1, false);
      Stk->noise = (double *)malloc(Stk->nl*sizeof(double));

    }else{
      if(Mpi->rank>0){  
        Stk->Jacobian = (double ***)TENSOR_DBL(0, 3, 1, \
            9, 0, Stk->nl-1, false); 
        Stk->noise = (double *)malloc(Stk->nl*sizeof(double));
      }else{
        free(Stk->Lambda);
      }
    }   

    Input->delta_v = (Stk->Lambda[Stk->nl-1]-Stk->Lambda[0])/Stk->nl \
        /Input->Lines[0].Lambda0*Par_C/1e3;

    if(Mpi->rank==0){
      sprintf(MeSS, "\n -- reading the wavelength -- \n");
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }

    return 0;

ERR1: 
    status = -1;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    Error(enum_error, rname, MeSS, NULL);
    return 0;

}

/*--------------------------------------------------------------------------------*/

extern int rprofile(STRUCT_INPUT *Input, int coord[], STRUCT_STK *Stk){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Save the best parameters.
      Record of revisions:
        17 Jun. 2024 (Hao Li)
      Input parameters:
        Input, the input configuration.
        coord, coordinates of the pixel.
        Stk, structure with Stokes profiles.
      Output parameters:
        Stk, structure with Stokes profiles.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    // CFITSIO status value MUST be initialized to zero!
    int status = 0, i;

    // get the first pixel
    Input->fpixel[2] = coord[0]+1;
    Input->fpixel[3] = coord[1]+1;
        
    // open the file.    
    fits_open_file(&(Input->fptr_data), Input->Data_Path, READONLY, &status);

    // read the profle
    if(Input->type_data == enum_dbl){
      fits_read_pix(Input->fptr_data, TDOUBLE, Input->fpixel, Stk->nl*4, \
          NULL, Stk->prof[0], NULL, &status);
    }else if(Input->type_data == enum_flt){
      fits_read_pix(Input->fptr_data, TFLOAT, Input->fpixel, Stk->nl*4, \
          NULL, Input->data_flt, NULL, &status);  
      for(i=0;i<Stk->nl*4;i++) Stk->prof[0][i] = (double)Input->data_flt[i];
     
    }else{
      fits_read_pix(Input->fptr_data, TSHORT, Input->fpixel, Stk->nl*4, \
          NULL, Input->data_int, NULL, &status);
      for(i=0;i<Stk->nl*4;i++) Stk->prof[0][i] = (double)Input->data_int[i];
    }

////////////////////SDO filter
/*
    for(i=0;i<4;i++){
      Stk->prof[i][0] *= 1.49279/1.43858806880917;
      Stk->prof[i][1] *= 1.49279/1.49139978559615;
      Stk->prof[i][2] *= 1.49279/1.52324167542270;
      Stk->prof[i][3] *= 1.49279/1.52811487224149;
      Stk->prof[i][4] *= 1.49279/1.50984871281028;
      Stk->prof[i][5] *= 1.49279/1.46553486521323;
    }
*/
    // close the file.
    fits_close_file(Input->fptr_data, &status);

    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int rprofileall(STRUCT_INPUT *Input, STRUCT_STK *Stk){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Save the best parameters.
      Record of revisions:
        27 Nov. 2024 (Hao Li)
      Input parameters:
        Input, the input configuration.
        Stk, structure with Stokes profiles.
      Output parameters:
        Stk, structure with Stokes profiles.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    const char *rname = "rprofileall";

    // CFITSIO status value MUST be initialized to zero!
    int status = 0;

    // open the file.    
    fits_open_file(&(Input->fptr_data), Input->Data_Path, READONLY, &status);
    
    // the first pixel
    long fpixel[4];

    fpixel[0] = 1;
    fpixel[1] = 1;
    fpixel[2] = 1;
    fpixel[3] = 1;

    // read the profle
    if(Input->type_data == enum_flt){
      fits_read_pix(Input->fptr_data, TFLOAT, fpixel, \
          Stk->nl*4*Input->nx*Input->ny, NULL, Stk->profall[0][0], \
          NULL, &status);  
    }else{
      sprintf(MeSS, "only float is supported in the fast mode.\n");
      Error(enum_error, rname, MeSS, NULL);
      return 1;  
    }

////////////////////SDO filter
/*
    for(i=0;i<4;i++){
      Stk->prof[i][0] *= 1.49279/1.43858806880917;
      Stk->prof[i][1] *= 1.49279/1.49139978559615;
      Stk->prof[i][2] *= 1.49279/1.52324167542270;
      Stk->prof[i][3] *= 1.49279/1.52811487224149;
      Stk->prof[i][4] *= 1.49279/1.50984871281028;
      Stk->prof[i][5] *= 1.49279/1.46553486521323;
    }
*/
    // close the file.
    fits_close_file(Input->fptr_data, &status);

    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int **cache_init(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, \
        STRUCT_STK *Stk, int *status){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        initialize or read the cache and result files.
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        Input, the input configuration.
        Mpi, structure with mpi configuration.
        Stk, structure with Stokes profiles.
      Output parameters:
        status, status of the fits file mutipulation.
      Return:
        return the pointer to the cache.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    const char *rname = "cache_init";

    // CFITSIO status value MUST be initialized to zero!
    *status = 0;

    long naxes[4];
    bool lcache = false;
    int nx, ny, npar = 10, naxis, **cache = NULL;
    Input->nxcache = 1+Input->sol_box[0][1]-Input->sol_box[0][0];
    Input->nycache = 1+Input->sol_box[1][1]-Input->sol_box[1][0];

    Input->counts = Input->nxcache*Input->nycache;

    char filename[Max_Line_Length+1];

    if(Mpi->rank==0){

      cache = (int **)MATRIX(0, Input->nxcache-1, 0, Input->nycache-1, \
          enum_int, true);

      Input->cache = fopen(Input->Cache_Path, "rb");
      if(Input->cache != NULL){
        char tmp[5];
        fread(tmp,sizeof(char),4,Input->cache);
        tmp[4] = '\0';

        if(strcmp(tmp, "cach") != 0){
          sprintf(MeSS, " ** %.*s if not a cache file.\n", \
              (int)strlen(Input->Cache_Path), Input->Cache_Path);
          VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>1);
          lcache = false;
        }else{

          fread(&nx, sizeof(int), 1, Input->cache);
          fread(&ny, sizeof(int), 1, Input->cache);

          if(nx==Input->nxcache && ny==Input->nycache){

            fread(cache[0], sizeof(int), Input->nxcache*Input->nycache, \
                Input->cache);

            sprintf(MeSS, "\n *** cache file readed. ***\n");
            VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>1);
            lcache = true;

          }else{
            sprintf(MeSS, " ** error in cache file size.\n");
            VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>1);
            lcache = false;

          }
        }
        fclose(Input->cache);
      }

      if(!lcache){

        sprintf(MeSS, "\n  *** no correct cache file found. ***");
        VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>1);
        sprintf(MeSS, "\n   ** creating a cache file.");
        VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>1);

        Input->cache = fopen(Input->Cache_Path, "wb");
        char tmp[5] = "cache";
        fwrite(tmp,sizeof(char), 4, Input->cache);
      
        fwrite(&Input->nxcache, sizeof(int), 1, Input->cache);
        fwrite(&Input->nycache, sizeof(int), 1, Input->cache);

        fwrite(cache[0],sizeof(int), Input->nxcache*Input->nycache, \
            Input->cache);

        fclose(Input->cache);

        // create a result file
        sprintf(MeSS, "\n   ** creating a fits file for the result.");
        VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>1);


        sprintf(filename, "!%s",Input->Result_Path);
        fits_create_file(&(Input->fptr_res), filename, status);

        //fits_create_file(&(Input->fptr_res), Input->Result_Path, status);

        if(*status!=0){
          sprintf(MeSS, "error in creaing the result file.\n");
          goto ERR1;
        }

        naxis = 3;
        naxes[0] = npar;
        naxes[1] = Input->nxcache;
        naxes[2] = Input->nycache;

        fits_create_img(Input->fptr_res, DOUBLE_IMG, naxis, naxes, status);
        
        if(*status!=0){
          sprintf(MeSS, "error in creaing the img hdu.\n");
          goto ERR1;
        }

        fits_close_file(Input->fptr_res, status);

        if(*status!=0){
          sprintf(MeSS, "error in closing the result file.\n");
          goto ERR1;
        }

        if(!Input->output_fit){
          MPI_Bcast(status, 1, MPI_INT, 0, MPI_COMM_WORLD);
          return cache;
        }

        sprintf(MeSS, "\n   ** creating a fits file for " \
            "the fitting profile.");
        VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>1);


        sprintf(filename, "!%s",Input->Fit_Path);
        fits_create_file(&(Input->fptr_fit), filename, status);

        //fits_create_file(&(Input->fptr_fit), Input->Fit_Path, status);

        if(*status!=0){
          sprintf(MeSS, "error in creaing the fit file.\n");
          goto ERR1;
        }

        naxis = 4;
        naxes[0] = Stk->nl;
        naxes[1] = 4;
        naxes[2] = Input->nxcache;
        naxes[3] = Input->nycache;

        fits_create_img(Input->fptr_fit, DOUBLE_IMG, naxis, naxes, \
            status);

        if(*status!=0){
          sprintf(MeSS, "error in creaing the img hdu "\
              "of the fit file.\n");
          goto ERR1;
        }
        fits_close_file(Input->fptr_fit, status);

        if(*status!=0){
          sprintf(MeSS, "error in closing the fit file.\n");
          goto ERR1;
        }

      }else{

        sprintf(MeSS, "\n   ** trying to open the corresponding "\
            "result file.\n");
        VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>1);

        fits_open_file(&(Input->fptr_res), Input->Result_Path, \
            READONLY, status);

        if(*status!=0){
          sprintf(MeSS, "error in opening the result file.\n");
          goto ERR1;
        }

        // get the dimension of the data hdu
        fits_get_img_dim(Input->fptr_res, &naxis, status);

        if(*status!=0){
          sprintf(MeSS, "error in getting dimension of the"\
              " result file.\n");
          goto ERR1;
        }

        if(naxis!=3){ 
          sprintf(MeSS, "error in the dimension of the"\
              " result file.\n");
          goto ERR1;
        }

        // get the size of each dimension 
        fits_get_img_size(Input->fptr_res, naxis, naxes, status);

        if(*status!=0){
          sprintf(MeSS, "error in getting size of the"\
              " result file.\n");
          goto ERR1;
        }

        if(naxes[0]!=npar||naxes[1]!=Input->nxcache \
            ||naxes[2]!=Input->nycache){
          sprintf(MeSS, "error in the size of the"\
              " result file.\n");
          goto ERR1;
        } 
          
        fits_close_file(Input->fptr_res, status);

        if(*status!=0){
          sprintf(MeSS, "error in closing the result file.\n");
          goto ERR1;
        }

        if(!Input->output_fit) {
          MPI_Bcast(status, 1, MPI_INT, 0, MPI_COMM_WORLD);
          return cache;
        }

        fits_open_file(&(Input->fptr_fit), Input->Fit_Path, READONLY, \
            status);

        if(*status!=0){
          sprintf(MeSS, "error in opening the fit file.\n");
          goto ERR1;
        }

        // get the dimension of the data hdu
        fits_get_img_dim(Input->fptr_fit, &naxis, status);

        if(*status!=0){
          sprintf(MeSS, "error in getting dimension of the fit file.\n");
          goto ERR1;
        }

        if(naxis!=4){ 
          sprintf(MeSS, "error in the dimension of the fit file.\n");
          goto ERR1;         
        }

        fits_get_img_size(Input->fptr_fit, naxis, naxes, status);

        if(*status!=0){
          sprintf(MeSS, "error in getting size of the"\
              " fit file.\n");
          goto ERR1;
        }
 
        if(naxes[0]!=Stk->nl || naxes[1]!=4 || \
            naxes[2]!=Input->nxcache|| naxes[3] != Input->nycache){
          sprintf(MeSS, "error in the size of the fit file.\n");
          goto ERR1;
        }

        fits_close_file(Input->fptr_fit, status);

        if(*status!=0){
          sprintf(MeSS, "error in closing the fit file.\n");
          goto ERR1;
        }
  
      }

      MPI_Bcast(status, 1, MPI_INT, 0, MPI_COMM_WORLD);

    }else{

      MPI_Bcast(status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(*status!=0) ABORTED();
    
    }

    return cache;

ERR1: 
    *status = -1;
    FREE_MATRIX(cache, 0, 0, enum_int);
    MPI_Bcast(status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    Error(enum_error, rname, MeSS, NULL);

    return cache;

}

/*--------------------------------------------------------------------------------*/

extern int cache_write(STRUCT_INPUT *Input, int *coord){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Write the cache.
      Record of revisions:
        25 Apr. 2024 (Hao Li)
      Input parameters:
        Input, the input configuration.
        coord, coordinates of the pixel.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    int offset = 12+(coord[0]*Input->nycache+coord[1])*4;
    int one = 1;

    Input->cache = fopen("./cache", "r+b");

    fseek(Input->cache, offset, SEEK_SET);
    fwrite(&one, sizeof(int), 1, Input->cache);

    fclose(Input->cache);

    return 0;

}

/*--------------------------------------------------------------------------------*/
