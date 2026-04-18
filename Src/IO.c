
#include "IO.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        16 Apr. 2026  (Hao Li)
          --- Updates:  
              rename the keyword fastmode t cache_prof, and add an 
              addition keyword cache_inv to cache the inversion progress.

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

#define IO_VERBOSE(lvl, fmt, ...)                                      \
    do{                                                                \
      if(Input->verboselv>=lvl){                                       \
        snprintf(MeSS, sizeof(MeSS), fmt, ##__VA_ARGS__);              \
        LOG_WRITE(MeSS, true, true);                                   \
      }                                                                \
    }while(0)     

#ifdef USE_MPI                                                         
#define IO_ERROR(fmt, ...)                                             \
    do{                                                                \
      CLOSE_FILES();                                                   \
      int stat = -1;                                                     \
      MPI_Bcast(&stat, 1, MPI_INT, 0, MPI_COMM_WORLD);               \
      snprintf(MeSS, sizeof(MeSS), fmt, ##__VA_ARGS__);                \
      LOG_ERROR(ERR_LVL_ERROR, "IO", MeSS);                            \
    }while(0)
#else
#define IO_ERROR(fmt, ...)                                             \
    do{                                                                \
      CLOSE_FILES();                                                   \
      snprintf(MeSS, sizeof(MeSS), fmt, ##__VA_ARGS__);                \
      LOG_ERROR(ERR_LVL_ERROR, "IO", MeSS);                            \
    }while(0)
#endif  

static FILE *cache_file = NULL;
static fitsfile *obs_file = NULL;
static fitsfile *res_file = NULL;
static fitsfile *fit_file = NULL;
static long Fits_fpixel[4];

/*--------------------------------------------------------------------------------*/

static int CLOSE_FITS_FILE(fitsfile **fileptr){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        close a .fits file.
      Record of revisions:
        10 Mar. 2026.
      Input parameters:
        .
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(*fileptr){
      int status = 0;
      fits_close_file(*fileptr, &status);
      *fileptr = NULL;
      if(status){ 
        LOG_ERROR(ERR_LVL_WARNING, "CLOSE_FITS_FILE", \
          "Error closing FITS file \n");
        return status;
      }
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

static void CLOSE_STD_FILE(FILE **fileptr){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        close a file.
      Record of revisions:
        10 Mar. 2026.
      Input parameters:
        .
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(*fileptr){
      fclose(*fileptr);
      fileptr = NULL;
    }
}

/*--------------------------------------------------------------------------------*/

int rWavelength(STRUCT_INPUT *Input, STRUCT_STK *Stk, \
    STRUCT_MPI *Mpi){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        read the wavelength.
      Record of revisions:
        10 Mar. 2026.
      Input parameters:
        Input, the input configuration.
        Stk, structure with Stokes profiles.
        Mpi, structure with mpi configuration.
      Output parameters:
        Stk, structure with Stokes profiles.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    // CFITSIO status value MUST be initialized to zero!
    // number of the HDUs, the bit of each pixel, the number of axes
    int status = 0, bitpix, naxis;;
    // the length of each axis, the first pixel
    long naxes[4];

    Fits_fpixel[0] = 1;
    Fits_fpixel[1] = 1;

    fitsfile *fptr_wav;
  
    if(Mpi->rank==0){

      IO_VERBOSE(2, "\n -- reading the wavelength -- \n");

      fits_open_file(&obs_file, Input->Data_Path, \
          READONLY, &status);
      if(status) IO_ERROR("error in opening the data file: " \
          "status = %d \n", status);   

      fits_open_file(&fptr_wav, Input->Wav_Path, READONLY, &status);
      if(status) IO_ERROR("error in opening the wavelength file: " \
          "status = %d \n", status);   

      // get the dimension of the data hdu
      fits_get_img_dim(obs_file, &(Input->naxis), &status);
      if(status) IO_ERROR("error in getting the dimension of "\
          "the data hdu: status = %d \n", status);   

      if(Input->naxis<3 || Input->naxis>4){ 
        IO_ERROR("wrong dimension of the data hdu.\n"); 
      }

      // get the size of each dimension 
      fits_get_img_size(obs_file, Input->naxis, naxes, &status);
      if(status) IO_ERROR("error in getting the size of the data hdu: " \
          "status = %d \n", status);   
      if(naxes[1]!=4) IO_ERROR("wrong size of the data hdu.\n");         

      Stk->nw = naxes[0];
      Input->nx = naxes[2];
      if(Input->naxis==3){
        Input->ny = 1;        
      }else{
        Input->ny = naxes[3];
      }
    }

#ifdef USE_MPI
     MPI_Bcast(&Input->nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast(&Input->ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    // check the solution box
    if(Input->sol_box[0][0]<0||Input->sol_box[0][0]>=Input->nx) \
        Input->sol_box[0][0] = 0;
    if(Input->sol_box[0][1]<0||Input->sol_box[0][1]>=Input->nx) \
        Input->sol_box[0][1] = Input->nx-1;
    if(Input->sol_box[1][0]<0||Input->sol_box[1][0]>=Input->ny) \
        Input->sol_box[1][0] = 0;
    if(Input->sol_box[1][1]<0||Input->sol_box[1][1]>=Input->ny) \
        Input->sol_box[1][1] = Input->ny-1;

    if(Mpi->rank==0){

      // get the type of the hdu
      fits_get_img_type(obs_file, &bitpix, &status);
      if(status) IO_ERROR("error in getting the data type " \
          "of data hdu: status = %d \n", status);   
      
      if(bitpix==16){ 
        Input->type_data = enum_int;  
      }else if(bitpix==-32){
        Input->type_data = enum_flt;
      }else if(bitpix==-64){
        Input->type_data = enum_dbl;
      }else{
        IO_ERROR("error in the data type of data hdu: " \
            "bitpix = %d.\n", bitpix);
      }
      
      // get the dimension of the wavlength hdu
      fits_get_img_dim(fptr_wav, &naxis, &status);
      if(status) IO_ERROR("error in getting the dimension " \
          "of wavelength hdu: status = %d \n", status);   
      if(naxis!=1) IO_ERROR("error in the dimension of " \
          "wavelength hdu.\n");

      // get the size of the wavelength 
      fits_get_img_size(fptr_wav, naxis, naxes, &status);
      if(status) IO_ERROR("error in getting the size of " \
          "the wavelength: status = %d \n", status);   
      if(Stk->nw!=naxes[0]) IO_ERROR("wavelength is not " \
          "consistent with the profile.\n");
      
      // get the type of the hdu
      fits_get_img_type(fptr_wav, &bitpix, &status);
      if(status) IO_ERROR("error in getting the data type of " \
          "the wavelength hdu: status = %d \n", status);   

      Stk->Lambda = (double *)malloc(Stk->nw*sizeof(double));

      if(bitpix==-32 || bitpix==-64){
        fits_read_pix(fptr_wav, TDOUBLE, Fits_fpixel, Stk->nw, NULL, \
            Stk->Lambda, NULL, &status);     
      }else{
        free(Stk->Lambda);
        IO_ERROR("error in the data type of the wavelength hdu.\n");
      }
      if(status) IO_ERROR("error in reading the wavelength: " \
          "status = %d \n", status);   

      status = CLOSE_FITS_FILE(&fptr_wav);
      fptr_wav = NULL;
      if(status) IO_ERROR("error in closing the wavelength file: " \
          "status = %d \n", status);   

      status = CLOSE_FITS_FILE(&obs_file);
      if(status) IO_ERROR("error in closing the data file: " \
          "status = %d \n", status);

#ifdef USE_MPI
      // broadcast to slaves
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&Stk->nw, 1, MPI_INT, 0, MPI_COMM_WORLD);
      // broadcast the wavelengths
      MPI_Bcast(Stk->Lambda, Stk->nw, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }else{

      // if error in reading wavelength
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(status) ABORTED();

      //  receive the size of the wavelength
      MPI_Bcast(&Stk->nw, 1, MPI_INT, 0, MPI_COMM_WORLD);
      Stk->Lambda = (double *)malloc(Stk->nw*sizeof(double));

      // receive the wavelength
      MPI_Bcast(Stk->Lambda, Stk->nw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    }

    

    return 0;
}

/*--------------------------------------------------------------------------------*/

int rProfile(STRUCT_INPUT *Input, STRUCT_STK *Stk, STRUCT_SUBSET *Subset){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        read profiles.
      Record of revisions:
        10 Mar. 2026.
      Input parameters:
        Input, the input configuration.
        coord, coordinates of the pixel.
        Stk, structure with Stokes profiles.
        Subset, a structure storing the pixels to read.
      Output parameters:
        Stk, structure with Stokes profiles.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    int status = 0;

    int nx = Input->cache_header.nx;

    // get the first pixel
    Fits_fpixel[2] = Subset->coord[0]+1;
    Fits_fpixel[3] = Subset->coord[1]+1;

    if(Subset->nProf >= Input->counts-Subset->pcounts){
      Subset->nProf = Input->counts-Subset->pcounts;
    }
    Subset->pcounts += Subset->nProf;
      
    int delta1 = Input->sol_box[0][1]-Subset->coord[0]+1;
    int delta = 0, nline = 0, delta2 = 0;

    if(Subset->nProf<delta1){
      Subset->coord[0] += Subset->nProf;
    }else if(Subset->nProf==delta1){

      Subset->coord[0] = Input->sol_box[0][0];
      Subset->coord[1]++;
    }else{
      delta = Subset->nProf-delta1;
      nline = delta/nx;
      delta2 = delta-nline*nx;
      Subset->coord[0] = Input->sol_box[0][0]+delta2;
      Subset->coord[1] += 1+nline;
    }

    fits_open_file(&obs_file, Input->Data_Path, \
        READONLY, &status);
    if(status) IO_ERROR("error in opening the data file: " \
        "status = %d \n", status); 

    if(Input->Subset_flg || Subset->nProf<=delta1){
      fits_read_pix(obs_file, TDOUBLE, Fits_fpixel, Subset->nProf*Stk->nw*4, \
          NULL, Input->profbuff, NULL, &status);

    }else{
      double *ptr = Input->profbuff;
   
      fits_read_pix(obs_file, TDOUBLE, Fits_fpixel, \
          delta1*Stk->nw*4, NULL, ptr, NULL, &status);        
      ptr += delta1*Stk->nw*4;

      Fits_fpixel[2] = Input->sol_box[0][0]+1;
      Fits_fpixel[3]++;

      for(int il=0; il<nline; il++){
        fits_read_pix(obs_file, TDOUBLE, Fits_fpixel, nx*Stk->nw*4, \
           NULL, ptr, NULL, &status);  
        ptr += nx*Stk->nw*4;   
        Fits_fpixel[3]++;
      }
      if(delta2>0){
        fits_read_pix(obs_file, TDOUBLE, Fits_fpixel, \
            delta2*Stk->nw*4, NULL, ptr, NULL, &status);  
        ptr += delta2*Stk->nw*4;   
        
      }
    }

    if(status) IO_ERROR("error in reading profiles: " \
        "status = %d \n", status);  

    status = CLOSE_FITS_FILE(&obs_file);
    if(status) IO_ERROR("rProfile error in closing the data file: " \
        "status = %d \n", status);  

    return status;
}

/*--------------------------------------------------------------------------------*/

int rProfilesll(STRUCT_INPUT *Input, STRUCT_STK *Stk){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        read all profiles.
      Record of revisions:
        16 Mar. 2026. 
      Input parameters:
        Input, the input configuration.
        Stk, structure with Stokes profiles.
      Output parameters:
        Stk, structure with Stokes profiles.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    int status = 0;
    long fpixel[4], lpixel[4], inc[4] = {1,1,1,1};

    fpixel[0] = 1;
    fpixel[1] = 1;
    fpixel[2] = Input->sol_box[0][0]+1;
    fpixel[3] = Input->sol_box[1][0]+1;

    lpixel[0] = Stk->nw;
    lpixel[1] = 4;
    lpixel[2] = Input->sol_box[0][1]+1;
    lpixel[3] = Input->sol_box[1][1]+1;

    fits_open_file(&obs_file, Input->Data_Path, \
        READONLY, &status);
    if(status) IO_ERROR("error in opening the data file: " \
        "status = %d \n", status); 

    fits_read_subset(obs_file, TDOUBLE, fpixel, lpixel, inc, \
        NULL, Input->profbuff, NULL, &status);

    if(status) IO_ERROR("error in reading the subset: " \
        "status = %d \n", status);   

    status = CLOSE_FITS_FILE(&obs_file);
    if(status) IO_ERROR("error in closing the data file: " \
        "status = %d \n", status);   
         
    return status;
}

/*--------------------------------------------------------------------------------*/

int PixelMV(STRUCT_INPUT *Input, STRUCT_SUBSET *Subset){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        moving to the next subset.
      Record of revisions:
        10 Mar. 2026.
      Input parameters:
        Input, the input configuration.
        Subset, a structure storing the pixels to read.
      Output parameters:
        Subset, a structure storing the pixels to read.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(Subset->nProf >= Input->counts-Subset->pcounts){
      Subset->nProf = Input->counts-Subset->pcounts;

      Subset->coord[0] = Input->sol_box[0][0];
      Subset->coord[1] = Input->sol_box[1][0]+1;

      Subset->pcounts = Input->counts;
      return 0;
    }else{
      int n1 = Input->sol_box[0][1]-Subset->coord[0]+1;

      int delta1 = 0, nline = 0, delta2 = 0;
      Subset->pcounts += Subset->nProf;

      if(Subset->nProf<n1){
        Subset->coord[0] += Subset->nProf;
      }else if(Subset->nProf==n1){
        Subset->coord[0] = Input->sol_box[0][0];
        Subset->coord[1]++;

      }else{
        delta1 = Subset->nProf-n1;
        nline = delta1/Input->cache_header.nx;
        delta2 = delta1-nline*Input->cache_header.nx;
        Subset->coord[0] = Input->sol_box[0][0]+delta2;
        Subset->coord[1] += 1+nline;
      }
    } 
        
    return 1;
}

/*--------------------------------------------------------------------------------*/

int CACHE_INIT(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, \
    STRUCT_STK *Stk){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        initialize or read the cache and result files.
      Record of revisions:
        10 Mar. 2026.
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

    memcpy(Input->cache_header.magic,"cach",4);
    // size of the cache matrix.
    int nxcache = Input->cache_header.nx;
    int nycache = Input->cache_header.ny;
    int ncache = Input->cache_header.ncache;

    Input->counts = nxcache*nycache;

    bool lcache = false;
    STRUCT_CACHE *header = &(Input->cache_header);
    char filename[Max_Line_Length+2];

    // CFITSIO status value MUST be initialized to zero!
    int status = 0, npar = 10, naxis;
    long naxes[4];

    if(Mpi->rank==0){
      Input->cache = (int *)calloc(ncache, sizeof(int));
      if(Input->cache_prof || !Input->cache_inv){
        lcache = false;
      }else{
        cache_file = fopen(Input->Cache_Path, "rb");

        if(cache_file != NULL){
          fread(header, sizeof(STRUCT_CACHE), 1, cache_file);
          lcache = (memcmp(header->magic, "cach", 4) == 0);

          if(lcache && nxcache == header->nx && nycache == header->ny \
              && ncache == header->ncache){
            
            fread(Input->cache, sizeof(int), ncache, cache_file);
            IO_VERBOSE(2, "\n *** cache file read. ***\n");  

          }else{
            header->nx = nxcache;
            header->ny = nycache;
            lcache = false;
            IO_VERBOSE(2, "\n ** %s is not a correct cache file.", \
                Input->Cache_Path);     
          }
          CLOSE_STD_FILE(&cache_file);
        }
      }

      if(!lcache){
        if(!Input->cache_prof && Input->cache_inv){
          IO_VERBOSE(2, "\n  *** no correct cache file found. ***");     
          IO_VERBOSE(2, "\n   ** creating a cache file.");  

          cache_file = fopen(Input->Cache_Path, "wb");
          fwrite(header, sizeof(STRUCT_CACHE), 1, cache_file);
          fwrite(Input->cache, sizeof(int), ncache, cache_file);
          CLOSE_STD_FILE(&cache_file);
        }

        IO_VERBOSE(2, "\n   ** creating a fits file for the result.");

        snprintf(filename, sizeof(filename), "!%s", \
            Input->Result_Path);
        fits_create_file(&res_file, filename, &status);
        if(status) IO_ERROR("error in creaing the result file: " \
            "status = %d \n", status);   

        naxis = 3;
        naxes[0] = npar;
        naxes[1] = nxcache;
        naxes[2] = nycache;

        fits_create_img(res_file, DOUBLE_IMG, naxis, naxes, &status);
        if(status) IO_ERROR("error in creaing the img hdu: " \
            "status = %d \n", status);    

        status = CLOSE_FITS_FILE(&res_file);
        if(status) IO_ERROR("error in closing the result file.\n");

        if(Input->output_fit){

          IO_VERBOSE(2, "\n   ** creating a fits file for " \
              "the fitting profile.");

          snprintf(filename, sizeof(filename), "!%s", Input->Fit_Path);
          fits_create_file(&(fit_file), filename, &status);
          if(status) IO_ERROR("error in creaing the fit file: " \
            "status = %d \n", status);   

          naxis = 4;
          naxes[0] = Stk->nw;
          naxes[1] = 4;
          naxes[2] = nxcache;
          naxes[3] = nycache;

          fits_create_img(fit_file, DOUBLE_IMG, naxis, naxes, &status);
          if(status) IO_ERROR("error in creaing the img hdu "\
                "of the fit file: status = %d \n", status);   

          status = CLOSE_FITS_FILE(&fit_file);
          if(status) IO_ERROR("error in closing the fit file.\n");
        }

      }else{

        IO_VERBOSE(2, "\n   ** opening the corresponding "\
            "result files.");
        
        fits_open_file(&res_file, Input->Result_Path, \
            READWRITE, &status);
 
        if(status) IO_ERROR("error in opening the result file: " \
            "status = %d \n",status);

        // get the dimension of the data hdu
        fits_get_img_dim(res_file, &naxis, &status);
        if(status) IO_ERROR("error in getting dimension of the"\
            " result file: status = %d \n", status);   
        if(naxis!=3) IO_ERROR( "error in the dimension of the"\
            " result file: status = %d \n", status);   
        // get the size of each dimension 
        fits_get_img_size(res_file, naxis, naxes, &status);
        if(status) IO_ERROR("error in getting size of the"\
            " result file: status = %d \n", status);   

        if(naxes[0]!=npar || naxes[1]!=nxcache || naxes[2]!=nycache){
          IO_ERROR("wrong size of the result file: status = %d \n", \
              status);   
        } 
        status = CLOSE_FITS_FILE(&res_file);
        if(status) IO_ERROR("error in closing the result file.\n");

        if(Input->output_fit){

          fits_open_file(&fit_file, Input->Fit_Path, READWRITE, \
              &status);
          if(status) IO_ERROR("error in opening the fit file: " \
              "status = %d \n", status);   

          // get the dimension of the data hdu
          fits_get_img_dim(fit_file, &naxis, &status);
          if(status) IO_ERROR("error in getting dimension of" \
              " the fit file: status = %d \n", status);   
          if(naxis!=4) IO_ERROR( "error in the dimension " \
              "of the fit file.");   

          fits_get_img_size(fit_file, naxis, naxes, &status);
          if(status) IO_ERROR("error in getting size of the" \
              " fit file: status = %d \n", status);   
          
          if(naxes[0]!=Stk->nw || naxes[1]!=4 || \
              naxes[2]!=nxcache|| naxes[3] != nycache){
            IO_ERROR("wrong size of the fit file.\n");
          }

          status = CLOSE_FITS_FILE(&fit_file);
          if(status) IO_ERROR("error in closing the fit file.\n");

        }
      }

#ifdef USE_MPI
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);

    }else{

      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(status) ABORTED();
#endif 

    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

int WRITE_RESULT(STRUCT_INPUT *Input, STRUCT_SUBSET *Subset, STRUCT_STK *Stk){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Write the cache.
      Record of revisions:
        10 Mar. 2026.
      Input parameters:
        Input, the input configuration.
        Subset, a structure storing the pixels of reading.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    int status = 0;
    long fpixel[4];

    fpixel[0] = 1;
    fpixel[1] = Subset->coord[0]+1-Input->sol_box[0][0];
    fpixel[2] = Subset->coord[1]+1-Input->sol_box[1][0];

    fits_open_file(&res_file, Input->Result_Path, \
        READWRITE, &status);
    if(status) IO_ERROR("error in opening the result file: " \
        "status = %d \n",status);
    fits_write_pix(res_file, TDOUBLE, fpixel, 10*Subset->nProf, \
          Input->resbuff, &status);
    if(status) IO_ERROR("error in writing the result file: " \
        "status = %d \n", status);  

    status = CLOSE_FITS_FILE(&res_file);
    if(status) IO_ERROR("error in closing the result file.\n");

    if(Input->output_fit){
      fits_open_file(&fit_file, Input->Fit_Path, READWRITE, \
          &status);
      if(status) IO_ERROR("error in opening the fit file: " \
          "status = %d \n", status);   


      fpixel[1] = 1;
      fpixel[2] = Subset->coord[0]+1-Input->sol_box[0][0];
      fpixel[3] = Subset->coord[1]+1-Input->sol_box[1][0];

      fits_write_pix(fit_file, TDOUBLE, fpixel, 4*Stk->nw*Subset->nProf, \
        Input->fitbuff, &status);

      if(status) IO_ERROR("error in writing the fit file: " \
        "status = %d \n", status);  

      status = CLOSE_FITS_FILE(&fit_file);
      if(status) IO_ERROR("error in closing the fit file.\n");
    }
 

    if(!Input->cache_prof && Input->cache_inv){
      int offset = sizeof(Input->cache_header) \
          +((Subset->coord[1]-Input->sol_box[1][0])*Input->cache_header.nx \
          +(Subset->coord[0]-Input->sol_box[0][0]))/Subset->nProf*sizeof(int);
      int one = 1;

      cache_file = fopen(Input->Cache_Path, "r+b");
      fseek(cache_file, offset, SEEK_SET);

      fwrite(&one, sizeof(int), 1, cache_file);
    
      CLOSE_STD_FILE(&cache_file);
    }
    return 0;
}

/*--------------------------------------------------------------------------------*/

int CLOSE_FILES(void){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        close the files.
      Record of revisions:
        10 Mar. 2026.
      Input parameters:
        .
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    LOG_FINALIZE(); 

    CLOSE_STD_FILE(&cache_file);

    int status = 0;
    if(CLOSE_FITS_FILE(&obs_file)!=0) status++;

    if(CLOSE_FITS_FILE(&res_file)!=0) status++;

    if(CLOSE_FITS_FILE(&fit_file)!=0) status++;

    return status;
}