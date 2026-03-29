
#include "RINPUT.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        09 Mar. 2026  (Hao Li)
          --- Updates:  
              removed the support of solver for multi lines to improve. 
              macros are used.

        18 Apr. 2025
          --- Updates: locate the line center 
                       update a few default parameters (Hao Li)

        11 Apr. 2025
          --- bugfix:  does not save the best fit (Hao Li)

        14 Jan. 2025
          --- bugfix:  fix the format for sprintf (Hao Li)

        27 Nov. 2024
          --- Updates:  fast mode of the inversion (Hao Li)

        28 Jun. 2024
          --- Initial commit (Hao Li)
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

#define IS_YES(str) (strcasecmp(str, "YES")==0)

/*--------------------------------------------------------------------------------*/

static int Get_Keys(STRUCT_KEYS Keywords[], STRUCT_INPUT *Input, \
    STRUCT_MPI *Mpi){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Convert the keywords to inversion configuration.
      Record of revisions:
        04 Mar. 2026.
      Input parameters:
        Keywords, structure with the input configuration.
      Output parameters:
        Input, structure with the input information.
        Stk, structure with the Stokes profiles.
        LM, structure with the Levenberg-Marquardt configuration.
        Par, structure with model parameters.
        Mpi, structure with Mpi information.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    const char *rname = "Get_Keys";

    // Pi Ratio of circumference to diameter          
    #define L_Pi 3.14159265358979323846

    char parameter[Key_Length];
    int nread, indx;

    #define INPUT_VERBOSE(indx, fmt, ...)                                 \
      do{                                                                 \
        if(Keywords[indx].set && Input->verboselv>0 && Mpi->rank==0){     \
          snprintf(MeSS, sizeof(MeSS), fmt, ##__VA_ARGS__);               \
          LOG_WRITE(MeSS, true, true);                                    \
        }                                                                 \
      }while(0)

    #define INPUT_FLAG(indx, var, desc)                                   \
      do{                                                                 \
        STR_TOUPPER(Keywords[indx].line);                                 \
        var = IS_YES(Keywords[indx].line);                                \
        INPUT_VERBOSE(indx, "\n %s: %s", desc, var ? "Yes" : "No");    \
      }while(0)


    #define PAR_BOUNDS(indx, par_indx, desc, unit)                          \
      do{                                                                   \
        nread = sscanf(Keywords[indx].line, "%lf, %lf",                     \
            &Input->Limits[par_indx][0], &Input->Limits[par_indx][1]);      \
        if(nread!=2 && Mpi->rank==0){                                       \
          snprintf(MeSS, sizeof(MeSS),                                      \
              "Error reading two doubles for key %s.\n", desc);             \
          LOG_ERROR(ERR_LVL_ERROR, rname, MeSS);                            \
        }                                                                   \
        INPUT_VERBOSE(indx, "\n Bounds on %s: from %e to %e %s ", desc,     \
            Input->Limits[par_indx][0], Input->Limits[par_indx][1], unit);  \
      }while(0)

    Input->verboselv = atoi(Keywords[KEY_VERBOSE].line);
    INPUT_VERBOSE(KEY_VERBOSE, "\n Verbose level: %d", Input->verboselv);

    if(Input->verboselv && Mpi->rank!=0) LOG_INIT(Input->Log_Path);

    Input->Lambda0 = (double *)malloc(Input->nline*sizeof(double));
    Input->Geff = (double *)malloc(Input->nline*sizeof(double));
    for(int ii=0;ii<Input->nline;ii++){
      nread = sscanf(Input->tmp_lines[ii],"%lf, %lf", \
          &(Input->Lambda0[ii]), &(Input->Geff[ii]));

      if(nread!=2 && Mpi->rank== 0){
        LOG_ERROR(ERR_LVL_ERROR, rname, \
            "Error in reading the spectral line info.\n");
        return 0;   
      }
      INPUT_VERBOSE(KEY_LINES, \
          "\n Line center: %e, Effective Lande factor: %e", \
          Input->Lambda0[ii], Input->Geff[ii]);
    }


    Input->niter = atoi(Keywords[KEY_NUM_ITER].line);
    INPUT_VERBOSE(KEY_NUM_ITER, "\n Max number of iterations: %d", \
        Input->niter);


    Input->nruns = atoi(Keywords[KEY_NUM_RUN].line); 
    INPUT_VERBOSE(KEY_NUM_RUN, "\n Max number of inversions: %d", \
        Input->nruns);

    
    Input->Convg_Criteria = atof(Keywords[KEY_CONVG_CRITERIA].line);
    INPUT_VERBOSE(KEY_CONVG_CRITERIA, "\n Convergence criteria: %e", \
        Input->Convg_Criteria);
    
    
    Input->Chisq_Criteria = atof(Keywords[KEY_CHISQ_CRITERIA].line);
    INPUT_VERBOSE(KEY_CHISQ_CRITERIA, \
        "\n Chi-square criteria for redoing the inversion : %e", \
        Input->Chisq_Criteria);
    

    STR_COPY(Input->Data_Path, Max_Line_Length, \
        Keywords[KEY_DATA_PATH].line, \
        strlen(Keywords[KEY_DATA_PATH].line), true);
    STR_TRIM(Input->Data_Path);
    INPUT_VERBOSE(KEY_DATA_PATH, "\n Data file path: %s", \
        Input->Data_Path);
    if(!FILE_EXIST(Input->Data_Path) && Mpi->rank==0){
      LOG_ERROR(ERR_LVL_ERROR, rname, "Data file doesn't exist. \n");  
    }


    STR_COPY(Input->Wav_Path, Max_Line_Length, \
        Keywords[KEY_WAV_PATH].line, \
        strlen(Keywords[KEY_WAV_PATH].line), true);
    STR_TRIM(Input->Wav_Path);
    INPUT_VERBOSE(KEY_WAV_PATH, "\n Wavelength file path: %s", \
        Input->Wav_Path);
    if(!FILE_EXIST(Input->Wav_Path) && Mpi->rank==0){
      LOG_ERROR(ERR_LVL_ERROR, rname, \
          "wavelength file doesn't exist. \n");  
    }


    nread = sscanf(Keywords[KEY_WEIGHTS].line,"%lf, %lf, %lf, %lf", \
        &(Input->Weights[0]), &(Input->Weights[1]), \
        &(Input->Weights[2]), &(Input->Weights[3]));
    if(nread!=4 && Mpi->rank==0){
      LOG_ERROR(ERR_LVL_ERROR, rname, "4 floats needed for weights. \n");  
    }
    INPUT_VERBOSE(KEY_WEIGHTS, "\n Weights: %f %f %f %f", \
        Input->Weights[0], Input->Weights[1], Input->Weights[2], \
        Input->Weights[3]);


    nread = sscanf(Keywords[KEY_SOL_BOX].line,"%d, %d, %d, %d", \
        &(Input->sol_box[0][0]), &(Input->sol_box[0][1]), \
        &(Input->sol_box[1][0]), &(Input->sol_box[1][1]));
    if(nread!=4 && Mpi->rank==0){
      LOG_ERROR(ERR_LVL_ERROR, rname, \
          "4 integers needed for the sol_box. \n");  
    }
    INPUT_VERBOSE(KEY_SOL_BOX, \
        "\n Solution region\n X: %d to %d \n Y: %d to %d", \
        Input->sol_box[0][0], Input->sol_box[0][1], \
        Input->sol_box[1][0], Input->sol_box[1][1]);


    INPUT_FLAG(KEY_OUTPUT_FIT, Input->output_fit, \
        "Output fitting profiles");


    STR_COPY(Input->Cache_Path, Max_Line_Length, \
        Keywords[KEY_CACHE_PATH].line, \
        strlen(Keywords[KEY_CACHE_PATH].line), true);
    STR_TRIM(Input->Cache_Path);
    INPUT_VERBOSE(KEY_CACHE_PATH, "\n Cache path: %s", \
        Input->Cache_Path);


    STR_COPY(Input->Result_Path, Max_Line_Length, \
        Keywords[KEY_RESULT_PATH].line, \
        strlen(Keywords[KEY_RESULT_PATH].line), true);
    STR_TRIM(Input->Result_Path);
    INPUT_VERBOSE(KEY_RESULT_PATH, "\n Result path: %s", \
        Input->Result_Path);


    STR_COPY(Input->Fit_Path, Max_Line_Length, \
        Keywords[KEY_FIT_PATH].line, \
        strlen(Keywords[KEY_FIT_PATH].line), true);
    STR_TRIM(Input->Fit_Path);
    INPUT_VERBOSE(KEY_FIT_PATH, "\n Fit path: %s", Input->Fit_Path);


    Input->Damp = atof(Keywords[KEY_LM_DAMP].line);
    Input->Lam_accept = atof(Keywords[KEY_LM_ACCEPT].line);
    Input->Lam_reject = atof(Keywords[KEY_LM_REJECT].line);
    nread = sscanf(Keywords[KEY_LM_LIMITS].line,"%lf, %lf", \
        &(Input->Lam_Lim[0]), &(Input->Lam_Lim[1]));

    if(Input->Lam_accept<1.5 || Input->Lam_reject <1.5){
      LOG_ERROR(ERR_LVL_ERROR, rname, \
          "Lam_accept and Lam_reject should be larger than 1.5. \n");  
    }
    if(nread!=2 && Mpi->rank==0){
      LOG_ERROR(ERR_LVL_ERROR, rname, \
          "2 floats needed for the lambda limits. \n");  
    }
    INPUT_VERBOSE(KEY_LM_DAMP, \
        "\n LM damping: %e, accept: %e, reject: %e, limits: %e-%e", \
        Input->Damp, Input->Lam_accept, Input->Lam_reject, \
        Input->Lam_Lim[0], Input->Lam_Lim[1]);


    PAR_BOUNDS(KEY_BMOD_LIMIT, 0, "Bmod", "[G]");
    PAR_BOUNDS(KEY_BTHETA_LIMIT, 1, "Btheta", "[Pi]");  
    Input->Limits[1][0] *= L_Pi; 
    Input->Limits[1][1] *= L_Pi;
    PAR_BOUNDS(KEY_BPHI_LIMIT, 2, "Bphi", "[Pi]");    
    Input->Limits[2][0] *= L_Pi; 
    Input->Limits[2][1] *= L_Pi;
    PAR_BOUNDS(KEY_VLOS_LIMIT, 3, "Vlos", "[km/s]");
    PAR_BOUNDS(KEY_DOPPLER_LIMIT, 4, "Doppler width", "[mA]");
    PAR_BOUNDS(KEY_DAMPING_LIMIT, 5, "Damping width", "[Dopp]");
    PAR_BOUNDS(KEY_ETA_LIMIT, 6, "Eta", "");
    PAR_BOUNDS(KEY_CONT_LIMIT, 7, "Continuum", "[DN]");
    PAR_BOUNDS(KEY_BETA_LIMIT, 8, "Beta", "");


    Input->Regul = false;
    for(int ii=0;ii<9;ii++){
      indx = KEY_FIXED_BMOD+ii;
      STR_TOUPPER(Keywords[indx].line);
      if(IS_YES(Keywords[indx].line)){
        Input->inv[ii] = false;
        INPUT_VERBOSE(indx, "\n Parameter %d is fixed. \n", ii);
      }else{
        Input->inv[ii] = true;
      }
      indx = KEY_REGUL_BMOD+ii;

      nread = sscanf(Keywords[indx].line,"%[^,], %lf, %lf", parameter, \
          &(Input->value_const[ii]), &(Input->Regul_weight[ii]));

      Input->regl[ii] = false;
      if(nread == 3){
        STR_TRIM(parameter);
        STR_TOUPPER(parameter);
        if(IS_YES(parameter)){
          Input->regl[ii] = true;
          Input->Regul = true;
          if(ii==1||ii==2){
            Input->value_const[ii] *= L_Pi;
          }
          INPUT_VERBOSE(indx, \
              "\n Parameter %d: regularization constant %.6e, " \
              "weight %.6e ", ii, Input->value_const[ii], \
              Input->Regul_weight[ii]);
        }
      }
    }


    Input->VCoeffi = atof(Keywords[KEY_VCOEFFI].line);
    INPUT_VERBOSE(KEY_VCOEFFI, "\n VCoeffi: %e. ", Input->VCoeffi);


    Input->LCoeffi = atof(Keywords[KEY_LCOEFFI].line);
    INPUT_VERBOSE(KEY_LCOEFFI, "\n LCoeffi: %e", Input->LCoeffi);


    Input->Icriteria = atof(Keywords[KEY_ICRITERIA].line);
    INPUT_VERBOSE(KEY_ICRITERIA, "\n Intensity threshold: %e", \
        Input->Icriteria);


    Input->Threshold = atof(Keywords[KEY_SVD_THRESHOLD].line);
    INPUT_VERBOSE(KEY_SVD_THRESHOLD, "\n SVD threshold: %e", \
        Input->Threshold);
 

    INPUT_FLAG(KEY_HMIREF, Input->HMI_REF, "HMI reference direction");


    INPUT_FLAG(KEY_FASTMODE, Input->fastmode, "Fast mode");
    if(Input->fastmode) Input->output_fit = false;

    Mpi->nthreads = 1;
#ifdef USE_OPENMP
    Mpi->nthreads = atoi(Keywords[KEY_NTHREADS].line);
#ifdef USE_MPI
    if(Mpi->size>1 && Mpi->rank==0){
      Mpi->nthreads = 1;
    }
#endif
    if(Keywords[KEY_NTHREADS].set && Input->verboselv>0){     
      snprintf(MeSS, sizeof(MeSS), "\n Rank %d using %d OpenMP threads", \
          Mpi->rank, Mpi->nthreads);                
      LOG_WRITE(MeSS, true, true);                                    
    }   
#endif

    Input->nfigures = atoi(Keywords[KEY_NFIGURES].line); 
    INPUT_VERBOSE(KEY_NFIGURES, "\n Voigt function accuracy: %d", \
        Input->nfigures);

    nread = sscanf(Keywords[KEY_STEP].line,                           \
        "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf",                \
        &(Input->step[0]), &(Input->step[1]), &(Input->step[2]),      \
        &(Input->step[3]), &(Input->step[4]), &(Input->step[5]),      \
        &(Input->step[6]), &(Input->step[7]), &(Input->step[8]));
    if(nread!=9 && Mpi->rank==0){
      LOG_ERROR(ERR_LVL_ERROR, rname,                                 \
          "9 floats needed for the steps. \n");  
    }
    INPUT_VERBOSE(KEY_STEP,                                           \
        "\n LM steps: %e, %e, %e, %e, %e, %e, %e, %e, %e\n",          \
        Input->step[0], Input->step[1], Input->step[2],               \
        Input->step[3], Input->step[4], Input->step[5],               \
        Input->step[6], Input->step[7], Input->step[8]);


    Input->nProf = atoi(Keywords[KEY_NPROF].line); 
    INPUT_VERBOSE(KEY_NPROF, "\n profile number: %d", Input->nProf);

    return 0;
}

/*--------------------------------------------------------------------------------*/

int RDINPUT(const char Filename[], STRUCT_INPUT *Input, STRUCT_MPI *Mpi){
  
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Read the input.
      Record of revisions:
        04 Mar. 2026. 
      Input parameters:
        Filename, the path to the input file.
      Output parameters:
        Input, structure with the input information.
        Stk, structure with the Stokes profiles.
        LM, structure with the Levenberg-Marquardt configuration.
        Par, structure with model parameters.
        Mpi, structure with Mpi information.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    const char *rname = "RDINPUT";

    snprintf(Input->Log_Path, sizeof(Input->Log_Path), "./log_%05d", \
        Mpi->rank);

    if(Mpi->rank==0) LOG_INIT(Input->Log_Path);

    if(!FILE_EXIST(Filename) && Mpi->rank==0){
      snprintf(MeSS, sizeof(MeSS), "control file doesn't exist: %s \n", \
            Filename);
      LOG_ERROR(ERR_LVL_ERROR, rname, MeSS);
    }

    FILE *Fa = fopen(Filename, "r");
    if(!Fa && Mpi->rank==0){
      snprintf(MeSS, sizeof(MeSS), "Cannot open control file: %s \n", \
          Filename);
      LOG_ERROR(ERR_LVL_ERROR, rname, MeSS);
    }    

    STRUCT_KEYS Keywords[KEY_TOTAL] = {
      [KEY_LINES] = KEY_REQ("lines"),                                         //0
      [KEY_VERBOSE] = KEY_DEF("verbose","1"),                                 //1
      [KEY_NUM_ITER] = KEY_DEF("num_iter","60"),                              //2
      [KEY_NUM_RUN] = KEY_DEF("num_run","4"),                                 //3
      [KEY_CONVG_CRITERIA] = KEY_DEF("convg_criteria","1e-3"),                //4
      [KEY_CHISQ_CRITERIA] = KEY_DEF("chisq_criteria","1.0"),                 //5
      [KEY_DATA_PATH] = KEY_REQ("data_path"),                                 //6
      [KEY_WAV_PATH] = KEY_REQ("wavelength_path"),                            //7
      [KEY_WEIGHTS] =  KEY_DEF("weights","1., 1., 1., 1."),                   //8
      [KEY_SOL_BOX] = KEY_DEF("sol_box","-1, -1, -1, -1"),                    //9
      [KEY_OUTPUT_FIT] = KEY_DEF("output_fit","NO"),                          //10
      [KEY_CACHE_PATH] = KEY_DEF("cache_path","./cache"),                     //11
      [KEY_RESULT_PATH] = KEY_DEF("result_path","./res.fits"),                //12
      [KEY_FIT_PATH] = KEY_DEF("fit_path","./fit.fits"),                      //13
      [KEY_LM_DAMP] = KEY_DEF("damp","1e-3"),                                 //14
      [KEY_LM_ACCEPT] = KEY_DEF("lambda_accept","5.0"),                       //15
      [KEY_LM_REJECT] = KEY_DEF("lambda_reject","5.0"),                       //16
      [KEY_LM_LIMITS] = KEY_DEF("lambda_limits","1e-5, 1e3"),                 //17
      [KEY_BMOD_LIMIT] = KEY_DEF("Bounds_Bmod","5, 6000"),                    //18
      [KEY_BTHETA_LIMIT] = KEY_DEF("Bounds_Btheta","0, 1"),                   //19
      [KEY_BPHI_LIMIT] = KEY_DEF("Bounds_Bphi","0, 1"),                       //20
      [KEY_VLOS_LIMIT] = KEY_DEF("Bounds_Vlos","-20., 20."),                  //21
      [KEY_DOPPLER_LIMIT] = KEY_DEF("Bounds_Dopp","20, 70"),                  //22
      [KEY_DAMPING_LIMIT] = KEY_DEF("Bounds_Damp","0.4, 0.6"),                //23
      [KEY_ETA_LIMIT] = KEY_DEF("Bounds_Eta","2, 70"),                        //24
      [KEY_CONT_LIMIT] = KEY_DEF("Bounds_Continuum","0, 1e6"),                //25
      [KEY_BETA_LIMIT] = KEY_DEF("Bounds_Beta","0.1, 0.7"),                   //26
      [KEY_FIXED_BMOD] = KEY_DEF("Bmod_fixed","NO"),                          //27
      [KEY_FIXED_BTHETA] = KEY_DEF("Btheta_fixed","NO"),                      //28
      [KEY_FIXED_BPHI] = KEY_DEF("Bphi_fixed","NO"),                          //29
      [KEY_FIXED_VLOS] = KEY_DEF("Vlos_fixed","NO"),                          //30
      [KEY_FIXED_DOPPLER] = KEY_DEF("Dopp_fixed","NO"),                       //31
      [KEY_FIXED_DAMPING] = KEY_DEF("Damp_fixed","NO"),                       //32
      [KEY_FIXED_ETA] = KEY_DEF("Eta_fixed","NO"),                            //33
      [KEY_FIXED_CONT] = KEY_DEF("Continuum_fixed","NO"),                     //34
      [KEY_FIXED_BETA] = KEY_DEF("Beta_fixed","NO"),                          //35
      [KEY_REGUL_BMOD] = KEY_DEF("Bmod_Regularization","NO"),                 //36
      [KEY_REGUL_BTHETA] = KEY_DEF("Btheta_Regularization","YES, 0.5, 1e-4"), //37
      [KEY_REGUL_BPHI] = KEY_DEF("Bphi_Regularization","NO"),                 //38
      [KEY_REGUL_VLOS] = KEY_DEF("Vlos_Regularization","NO"),                 //39
      [KEY_REGUL_DOPPLER] = KEY_DEF("Dopp_Regularization","NO, 20, 0.05"),    //40
      [KEY_REGUL_DAMPING] = KEY_DEF("Damp_Regularization","NO, 0.5, 80"),     //41
      [KEY_REGUL_ETA] = KEY_DEF("Eta_Regularization","NO, 2.0, 0.01"),        //42
      [KEY_REGUL_CONT] = KEY_DEF("Continuum_Regularization","NO"),            //43
      [KEY_REGUL_BETA] =KEY_DEF("Beta_Regularization","NO"),                  //44
      [KEY_VCOEFFI] = KEY_DEF("VCoeffi","100."),                              //45
      [KEY_LCOEFFI] = KEY_DEF("LCoeffi","250."),                              //46
      [KEY_ICRITERIA] = KEY_DEF("Icriteria","100."),                          //47
      [KEY_SVD_THRESHOLD] = KEY_DEF("SVDthreshold","1e-6"),                   //48
      [KEY_HMIREF] = KEY_DEF("HMI_REF","No"),                                 //49
      [KEY_FASTMODE] = KEY_DEF("Fastmode","No"),                             //50
      [KEY_NTHREADS] = KEY_DEF("Nthreads","4"),                               //51
      [KEY_NFIGURES] = KEY_DEF("Nfigures","6"),                               //52
      [KEY_STEP] = KEY_DEF("Step", \
          "1000., 2., 2., 2., 20., 0.3, 30., 1e5, 0.3"),                      //53
      [KEY_NPROF] = KEY_DEF("nprofiles", "20")                                //54
    };

    char *lines = NULL, key[Key_Length], *value;
    size_t size = 0;
    int len_tot, len, nspace;
    Input->verboselv = 3;
    Input->nline = 0;
    bool found = false;
        
    while((len_tot=STR_READ_LINE(&lines, &size, Fa)) > 0){
      len = STR_INDEX_CHAR(lines, '=', 1);
      if(len <= 1) continue;

      if(len > Key_Length-1){
        if( Mpi->rank==0){
          snprintf(MeSS, sizeof(MeSS), \
              "\n keyword to long. \n %s \n line skipped. \n", \
              lines);
          LOG_ERROR(ERR_LVL_WARNING, rname, MeSS);
        }
        continue;
      }

      STR_COPY(key, Key_Length, lines, len, 1);
      STR_TRIM_RIGHT(key);
        
      value = lines+len+1;
      len_tot -= (len+1);
      found = false;
      nspace = STR_TRIM_LEFT(value);

      if(nspace>0) len_tot -= nspace;

      if(len_tot>Max_Line_Length-1){
        len_tot = Max_Line_Length-1;

        if(Mpi->rank==0){
          snprintf(MeSS, sizeof(MeSS), \
              "\n the value for keyword %s is too long: %s \n", \
              key, value);
          LOG_ERROR(ERR_LVL_WARNING, rname, MeSS);
        }
      }

      if(strcmp(key,Keywords[0].keyword)==0){
        STR_COPY(Input->tmp_lines[Input->nline], Max_Line_Length, value, \
            len_tot, false);
        Input->nline++;
        found = true;
        Keywords[0].set = true;

        if(Input->nline>NUM_LINES && Mpi->rank==0){
          snprintf(MeSS, sizeof(MeSS), \
              "\n keyword %s input more than %d \n", \
              key, NUM_LINES);
          LOG_ERROR(ERR_LVL_ERROR, rname, MeSS);
        }
      }else{
        for(int i=1; i<KEY_TOTAL; i++){
          if(strcmp(key,Keywords[i].keyword)==0){
            if(Keywords[i].set && Mpi->rank==0){
              snprintf(MeSS, sizeof(MeSS), \
                  "\n keyword %s redefined. Overwriting. \n", \
                  key);
              LOG_ERROR(ERR_LVL_WARNING, rname, MeSS);      
            }
            STR_COPY(Keywords[i].line, Max_Line_Length, value, len_tot, true);
            Keywords[i].set = true;

#ifdef DEBUG
            if(Mpi->rank==0){
              snprintf(MeSS, sizeof(MeSS), "\n   -- read keyword : %s - %s\n", \
                  Keywords[i].keyword, Keywords[i].line);
              LOG_WRITE(MeSS, false, true);
            }
#endif 

            found = true;
            break;
          }
        }
      }
      if(!found && Mpi->rank==0){
        snprintf(MeSS, sizeof(MeSS), "\n Unknown keyword: %s\n", key);
        LOG_ERROR(ERR_LVL_WARNING, rname, MeSS); 
      }
    }

    free(lines);
    fclose(Fa);

    for(int i=0; i<KEY_TOTAL; i++){
      if(Keywords[i].required && !Keywords[i].set && Mpi->rank==0){
        snprintf(MeSS, sizeof(MeSS), \
            "\n Required keyword missing: %.64s \n", Keywords[i].keyword);
        LOG_ERROR(ERR_LVL_ERROR, rname, MeSS);
      }
#ifdef DEBUG
      Keywords[i].set = true;
#endif
    }

    Get_Keys(Keywords, Input, Mpi);

    return 0;
}

/*--------------------------------------------------------------------------------*/
