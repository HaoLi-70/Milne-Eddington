
#include "READ_INPUT.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

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

static double Gfactor(double J,double L,double S);

static double Geffect(double Gu, double Gl, double Ju, double Jl);

static int Get_Keys(STRUCT_KEYS Keywords[], STRUCT_INPUT *Input, \
                    STRUCT_STK *Stk, STRUCT_LM *LM, \
                    STRUCT_PAR *Par, STRUCT_MPI *Mpi);

/*--------------------------------------------------------------------------------*/

static double Gfactor(double J,double L,double S){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Computes the Lande factor.
      Modified:
        25 Apr. 2018 (Hao Li)
      Input parameters:
        J, The total angular momentum of the electronic cloud.
        L, The total orbital angular momentum of the electronic cloud.
        S, The total spin of the electronic cloud.
      Return:
        Gfactor, The Lande factor.
      Method:
        L-S coupling asumption, and if the number is 0, return 0 for 
        convernience.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
    
    if(J==0) return 0;
        
    return 1.0+(J*(J+1.0)-L*(L+1.0)+S*(S+1.0))/(2.0*J*(J+1.0));
  
}

/*--------------------------------------------------------------------------------*/

static double Geffect(double Gu, double Gl, double Ju, double Jl){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Computes the effective Lande factor.
      Modified:
        27 Nov. 2024 (Hao Li)
      Input parameters:
        Gu, The lande factor of upper level.
        Gl, The lande factor of lower level.
        Ju, The total angular momentum of upper level.
        Jl, The total angular momentum of lower level.
      Return:
        Geffect, The effect Lande factor.
      References:
        LL04 Chapter 3, Equation 3.44.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
         
    return 0.5*(Gu+Gl)+0.25*(Gu-Gl)*(Ju*(Ju+1)-Jl*(Jl+1));
  
}

/*--------------------------------------------------------------------------------*/

static int Get_Keys(STRUCT_KEYS Keywords[], STRUCT_INPUT *Input, \
                    STRUCT_STK *Stk, STRUCT_LM *LM, \
                    STRUCT_PAR *Par, STRUCT_MPI *Mpi){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Convert the keywords to inversion configuration.
      Record of revisions:
        14 Jan. 2025 (Hao Li)
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
    char *token, parameter[Key_Length], *p;
    int indx, nread, len, i, itmp;
    long len_tot;
    

    indx = 1;
    Input->verboselv = atoi(Keywords[indx].line);
    if(Input->verboselv > 0){
      sprintf(MeSS, "\n verbose level: %d", Input->verboselv);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }
    

    indx = 0;
    Input->Lines = (STRUCT_MELINE *)malloc(Input->nline \
        *sizeof(STRUCT_MELINE));
    for(i=0;i<Input->nline;i++){
      nread = sscanf(Input->tmp_lines[i],"%lf, %lf", \
          &(Input->Lines[i].Lambda0), &(Input->Lines[i].Geff));

      if(nread!=2){
        sprintf(MeSS, "error in reading the line information.\n");
        Error(enum_error, rname, MeSS, NULL);
        return 0;   
      }else{
        if(Mpi->rank== 0){
          sprintf(MeSS, "\n line center: %e. effect Lande factor: %e", \
              Input->Lines[i].Lambda0, Input->Lines[i].Geff);
          VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
        }
      }
    }


    indx = 2;
    Input->niter = atoi(Keywords[indx].line);
    if(Input->verboselv > 0 && Mpi->rank == 0){
      sprintf(MeSS, "\n number of iteration : %d", Input->niter);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 3;
    Input->nrun = atoi(Keywords[indx].line);
    if(Input->verboselv > 0 && Mpi->rank == 0){
      sprintf(MeSS, "\n max number of inversion : %d", Input->nrun);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }
    
    
    indx = 4;
    Input->Convg_Criteria = atof(Keywords[indx].line);
    if(Input->verboselv > 0 && Mpi->rank == 0){
      sprintf(MeSS, "\n criteria for convergence: %e", \
          Input->Convg_Criteria);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);

    }
    
    
    indx = 5;
    Input->Chisq_Criteria = atof(Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n criteria for redoing the inversion: %e", \
          Input->Chisq_Criteria);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }
    

    indx = 6;
    String_Copy(Input->Data_Path, Keywords[indx].line, \
        strlen(Keywords[indx].line), true);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n path to the data file : %s", Input->Data_Path);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }
    if(!FILE_EXIST(Input->Data_Path)){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "data file doesn't exist.", \
            Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }


    indx = 7;
    String_Copy(Input->Wav_Path, Keywords[indx].line, \
        strlen(Keywords[indx].line), true);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n path to the wavelength file : %s", Input->Wav_Path);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }
    if(!FILE_EXIST(Input->Wav_Path)){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "wavelength file doesn't exist.", \
            Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }


    indx = 8;
    nread = sscanf(Keywords[indx].line,"%lf, %lf, %lf, %lf", \
                      &(Stk->Weights[0]), &(Stk->Weights[1]), \
                      &(Stk->Weights[2]), &(Stk->Weights[3]));
    if(nread!=4){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "4 floats needed for the weights.", \
            Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    Stk->Weights_SQR[0] = Stk->Weights[0]*Stk->Weights[0];
    Stk->Weights_SQR[1] = Stk->Weights[1]*Stk->Weights[1];
    Stk->Weights_SQR[2] = Stk->Weights[2]*Stk->Weights[2];
    Stk->Weights_SQR[3] = Stk->Weights[3]*Stk->Weights[3];
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n weights : %f %f %f %f\n", (Stk->Weights[0]), \
          (Stk->Weights[1]), (Stk->Weights[2]), (Stk->Weights[3]));
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 9;
    nread = sscanf(Keywords[indx].line,"%d, %d, %d, %d", \
                      &(Input->sol_box[0][0]), &(Input->sol_box[0][1]), \
                      &(Input->sol_box[1][0]), &(Input->sol_box[1][1]));
    if(nread!=4){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "4 integers needed for the sol_box.", \
            Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n solution region :");
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
      sprintf(MeSS, "\n     X from %d to %d", Input->sol_box[0][0], \
          Input->sol_box[0][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
      sprintf(MeSS, "\n     Y from %d to %d", Input->sol_box[0][0], \
          Input->sol_box[0][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 10;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->output_fit = true;
    }else{
      Input->output_fit = false;
    }
    if(Mpi->rank== 0){
      if(Input->output_fit){
        sprintf(MeSS, "\n output the fitting profiles : Yes");
        VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
      }
    }


    indx = 11;
    sprintf(Input->Cache_Path, "%s", Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n path to the cache file : %s", \
          Input->Cache_Path);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 12;
    sprintf(Input->Result_Path, "%s", Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n path to the result file : %s", \
          Input->Result_Path);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 13;
    sprintf(Input->Fit_Path, "%s", Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n path to the fit profile : %s", \
          Input->Fit_Path);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 14;
    LM->Damp = atof(Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n damping factor: %e", LM->Damp);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 15;
    LM->Lam_accept = atof(Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n factor to decrease the damping: %e", \
          LM->Lam_accept);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }
    
    
    indx = 16;
    LM->Lam_reject = atof(Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n factor to increase the dammping: %e", \
          LM->Lam_reject);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 17;
    nread = sscanf(Keywords[indx].line,"%lf, %lf", &(LM->Lam_Lim[0]), \
                      &(LM->Lam_Lim[1]));
    if(nread!=2){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "2 floats needed for the limits.", \
            Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n limits on the lambda factor: from %e to %e", \
          LM->Lam_Lim[0], LM->Lam_Lim[1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 18;
    nread = sscanf(Keywords[indx].line,"%lf, %lf", &(Par->Limits[1][0]), \
                      &(Par->Limits[1][1]));
    if(nread!=2){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "2 floats needed for the bounds of " \
            "the Bmod.", Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n bounds of the Bmod: from %e [G] to %e [G]", \
          Par->Limits[1][0], Par->Limits[1][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 19;
    nread = sscanf(Keywords[indx].line,"%lf, %lf", &(Par->Limits[2][0]), \
                      &(Par->Limits[2][1]));
    Par->Limits[2][0] *= Par_Pi;
    Par->Limits[2][1] *= Par_Pi;
    if(nread!=2){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "2 floats needed for the bounds of " \
            "the Btheta.", Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n bounds on the Btheta: from %e [Pi] to %e [Pi]", \
          Par->Limits[2][0], Par->Limits[2][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 20;
    nread = sscanf(Keywords[indx].line,"%lf, %lf", &(Par->Limits[3][0]), \
                      &(Par->Limits[3][1]));
    Par->Limits[3][0] *= Par_Pi;
    Par->Limits[3][1] *= Par_Pi;
    if(nread!=2){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "2 floats needed for the bounds of " \
            "the Bphi.", Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n bounds on the Bphi: from %e [Pi] to %e [Pi]", \
          Par->Limits[3][0], Par->Limits[3][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 21;
    nread = sscanf(Keywords[indx].line,"%lf, %lf", &(Par->Limits[4][0]), \
                      &(Par->Limits[4][1]));
    if(nread!=2){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "2 floats needed for the bounds of " \
            "the Vlos.", Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n bounds on the Vlos: from %e [km/s] to %e [km/s]", \
          Par->Limits[4][0], Par->Limits[4][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 22;
    nread = sscanf(Keywords[indx].line,"%lf, %lf", &(Par->Limits[5][0]), \
                      &(Par->Limits[5][1]));
    if(nread!=2){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "2 floats needed for the bounds of " \
            "the Doppler width.", Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n bounds on the Doppler width: from %e [mA] to %e [mA]", \
          Par->Limits[5][0], Par->Limits[5][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 23;
    nread = sscanf(Keywords[indx].line,"%lf, %lf", &(Par->Limits[6][0]), \
                      &(Par->Limits[6][1]));
    if(nread!=2){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "2 floats needed for the bounds of " \
            "the Damping width.", Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n bounds on the Damping width: from %e [Dopp] to %e [Dopp]", \
          Par->Limits[6][0], Par->Limits[6][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 24;
    nread = sscanf(Keywords[indx].line,"%lf, %lf", &(Par->Limits[7][0]), \
                      &(Par->Limits[7][1]));
    if(nread!=2){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "2 floats needed for the bounds of " \
            "the Eta.", Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n bounds on the Eta: from %e to %e", \
          Par->Limits[7][0], Par->Limits[7][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 25;
    nread = sscanf(Keywords[indx].line,"%lf, %lf", &(Par->Limits[8][0]), \
                      &(Par->Limits[8][1]));
    if(nread!=2){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "2 floats needed for the bounds of " \
            "the Continuum.", Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n bounds on the Continuum: from %e [DN] to %e [DN]", \
          Par->Limits[8][0], Par->Limits[8][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 26;
    nread = sscanf(Keywords[indx].line,"%lf, %lf", &(Par->Limits[9][0]), \
                      &(Par->Limits[9][1]));
    if(nread!=2){
      if (Mpi->rank == 0){
        Error(enum_error, rname, "2 floats needed for the bounds of " \
            "the Beta.", Input->Verbose_Path);
      }else{
        ABORTED();
      }
    }
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n bounds on the Beta: from %e to %e", \
          Par->Limits[9][0], Par->Limits[9][1]);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }
    

    itmp = indx;
    for(i=1;i<10;i++){
      indx = itmp+i;
      String_to_Upper(Keywords[indx].line);
      if(strcmp(Keywords[indx].line, "YES") == 0){
        Input->inv[i] = false;
      }else{
        Input->inv[i] = true;
      }
      if(Mpi->rank== 0){
        if(!Input->inv[i]){
          sprintf(MeSS, "\n Parameter %d is fixed. ", i);
          VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
        }
      }
    }

 
    itmp = 35;
    Input->Regul = false;
    for(i=1;i<10;i++){
      Input->regl[i] = false;
      indx = itmp+i;

      len_tot = strlen(Keywords[indx].line);
      len = Indx_Char(Keywords[indx].line, ',', 1);
      if(len > 1){
        len_tot -= len;
        
        String_Copy(parameter, Keywords[indx].line, len-1, 1);
        Trim(parameter, 3);
        String_to_Upper(parameter);
        if(strcmp(parameter, "YES") == 0){

          p = Keywords[indx].line+len;
          nread = sscanf(p,"%lf, %lf", &(Input->value_const[i]), \
            &(Input->Regul_weight[i]));

          if(nread==2){
            Input->regl[i] = true;
            Input->Regul = true;
            if(i==2||i==3){
              Input->value_const[i] *= Par_Pi;
            }
            if(Mpi->rank== 0){
              sprintf(MeSS, "\n Parameter %d: regularization constant %.6e; " \
                  " regularization weights: %.6e; ", i, Input->value_const[i], \
                  Input->Regul_weight[i]);
              VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);

            }
          }
        }
      }
    }
    

    indx = 45;
    Input->VCoeffi = atof(Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n coefficient for the initial guess of Blos : %e", \
          Input->VCoeffi);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 46;
    Input->LCoeffi = atof(Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n coefficient for the initial guess of Bpos : %e", \
          Input->LCoeffi);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 47;
    Input->Icriteria = atof(Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n Intensity threshold for the inversion : %e", \
          Input->Icriteria);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 48;
    Input->Threshold = atof(Keywords[indx].line);
    if(Mpi->rank== 0){
      sprintf(MeSS, "\n threshold for SVD : %e", \
          Input->Threshold);
      VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
    }


    indx = 49;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->Azimuth_Rotate = true;
    }else{
      Input->Azimuth_Rotate = false;
    }
    if(Mpi->rank== 0){
      if(Input->Azimuth_Rotate){
        sprintf(MeSS, "\n Rotate the azimuth : Yes");
        VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
      }
    }


    indx = 50;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->fastmode = true;
    }else{
      Input->fastmode = false;
    }
    if(Mpi->rank== 0){
      if(Input->Azimuth_Rotate){
        sprintf(MeSS, "\n Fast mode : Yes");
        VerboseM(MeSS, Input->Verbose_Path, Input->verboselv>0);
      }
    }

    return 0;

}

/*--------------------------------------------------------------------------------*/

extern int RDINPUT(char Filename[], STRUCT_INPUT *Input,  \
    STRUCT_STK *Stk, STRUCT_LM *LM, STRUCT_PAR *Par, STRUCT_MPI *Mpi){
  
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Read the input.
      Record of revisions:
        27 Nov. 2024 (Hao Li)  
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

    FILE *Fa = fopen(Filename, "r");
    sprintf(Input->Verbose_Path, "./verbose_%05d", Mpi->rank);
    int nKeywords;

    STRUCT_KEYS Keywords[] ={
      {"lines", "", false, true},  //0
      {"verbose", "1", false, false}, //1
      {"num_iter", "30", false, false}, //2
      {"num_run", "3", false, false}, //3
      {"convg_criteria", "1e-3", false, false}, //4
      {"chisq_criteria", "1.0", false, false}, //5
      {"data_path", "", false, true}, //6
      {"wavelength_path", "", false, true}, //7
      {"weights", "1., 1., 1., 1.", false, false},  //8
      {"sol_box", "-1, -1, -1, -1", false, false},  //9
      {"output_fit", "NO", false, false},  //10
      {"cache_path", "./cache", false, false}, //11
      {"result_path", "./res.fits", false, false}, //12
      {"fit_path", "./fit.fits", false, false}, //13
      {"damp", "1e-3", false, false}, //14
      {"lambda_accept", "5.0", false, false}, //15
      {"lambda_reject", "5.0", false, false}, //16
      {"lambda_limits", "1e-5, 1e3", false, false}, //17
      {"Bounds_Bmod", "5, 4000", false, false}, //18
      {"Bounds_Btheta", "0, 1", false, false}, //19
      {"Bounds_Bphi", "-2, 4", false, false}, //20
      {"Bounds_Vlos", "-20., 20.", false, false}, //21
      {"Bounds_Dopp", "20, 70", false, false}, //22
      {"Bounds_Damp", "0.4, 0.6", false, false}, //23
      {"Bounds_Eta", "2, 70", false, false}, //24
      {"Bounds_Continuum", "0, 1e6", false, false}, //25
      {"Bounds_Beta", "0.1, 0.7", false, false}, //26
      {"Bmod_fixed", "NO", false, false}, //27
      {"Btheta_fixed", "NO", false, false}, //28
      {"Bphi_fixed", "NO", false, false}, //29
      {"Vlos_fixed", "NO", false, false}, //30
      {"Dopp_fixed", "NO", false, false}, //31
      {"Damp_fixed", "NO", false, false}, //32
      {"Eta_fixed", "NO", false, false}, //33
      {"Continuum_fixed", "NO", false, false}, //34
      {"Beta_fixed", "NO", false, false}, //35
      {"Bmod_Regularization", "NO", false, false}, //36
      {"Btheta_Regularization", "YES, 0.5, 1e-4", false, false}, //37
      {"Bphi_Regularization", "NO", false, false}, //38
      {"Vlos_Regularization", "NO", false, false}, //39
      {"Dopp_Regularization", "NO, 20, 0.05", false, false}, //40
      {"Damp_Regularization", "NO, 0.5, 80", false, false}, //41
      {"Eta_Regularization", "NO, 2.0, 0.01", false, false}, //42
      {"Continuum_Regularization", "NO", false, false}, //43
      {"Beta_Regularization", "NO", false, false}, //44
      {"VCoeffi", "100.", false, false}, //45
      {"LCoeffi", "250.", false, false}, //46
      {"Icriteria", "1000.", false, false}, //47
      {"SVDthreshold", "1e-3", false, false}, //48
      {"Azm_Rot", "No", false, false}, //49
      {"Fastmode", "Yes", false, false} //50
    };
    nKeywords = sizeof(Keywords)/sizeof(STRUCT_KEYS);
    
    char lines[Max_Line_Length], key[Key_Length], *p;
    int len, i;
    long len_tot;
    bool neglect = true;

    Input->verboselv = 3;
    Input->nline = 0;

    while(Read_line(lines, Fa) > 0 ){
      len_tot = strlen(lines);
      len = Indx_Char(lines, '=', 1);
      if(len > 1){
        len_tot -= len;
        p = lines+len;
        String_Copy(key, lines, len-1, 1);
        Trim(key, 3);
        
        neglect = true;
        if(strcmp(key,Keywords[0].keyword)==0){
          String_Copy(Input->tmp_lines[Input->nline], p, len_tot, false);
          Input->nline++;
          neglect = false;
          Keywords[0].Set = true;
        }else{
          for(i=1; i<nKeywords; i++){

            if(strcmp(key,Keywords[i].keyword)==0){

              if(Keywords[i].Set == true && Mpi->rank == 0){
                sprintf(MeSS, "\n The keyword %.*s is set.", \
                    (int)(strlen(Keywords[i].keyword)), Keywords[i].keyword);
                Error(enum_warning, rname, MeSS, Input->Verbose_Path);
              }

              String_Copy(Keywords[i].line, p, len_tot, true);
              Trim(Keywords[i].line, 3);
              Keywords[i].Set = true;

              if(Mpi->rank == 0){
                sprintf(MeSS,"\n   -- read keyword : %.*s - %.*s", 
                    (int)(strlen(Keywords[i].keyword)), \
                    Keywords[i].keyword, (int)(strlen(Keywords[i].line)), \
                    Keywords[i].line);
                VerboseM(MeSS, Input->Verbose_Path, false);
              }

              neglect = false;
              break;
            }
          }
        }
        if(Mpi->rank == 0 && neglect){
          sprintf(MeSS, "\n %.*s is not keywords.", (int)(strlen(key)), key);
          Error(enum_warning, rname, MeSS, Input->Verbose_Path);
        }
      }
    }

    for(i=0; i<nKeywords; i++){
      if(Keywords[i].Required && !Keywords[i].Set){
        if(Mpi->rank == 0){
          sprintf(MeSS, "\n The keyword %.*s is required, but not set.\n", \
                  (int)(strlen(Keywords[i].keyword)), Keywords[i].keyword);
          Error(enum_error, rname, MeSS, Input->Verbose_Path);
        }else{
          ABORTED();
        }
      }
    }

    Par->Par = (double *)VECTOR(1, 9, enum_dbl, false);
    Par->Par_Best = (double *)VECTOR(1, 9, enum_dbl, false);
    Par->Par_tmp = (double *)VECTOR(1, 9, enum_dbl, false);

    Get_Keys(Keywords, Input, Stk, LM, Par, Mpi);
   
    /*
    double gu = Gfactor(0.0, 2.0, 2.0);
    double gl = Gfactor(1.0, 1.0, 2.0);
    Lines->Geffect = Geffect(gu, gl, 0.0, 1.0);
    */
    // Lines->Lambda in A, Par->Par[4] in mA, Lines->BShift in mA

    for(i=0;i<Input->nline;i++){
      Input->Lines[i].BShift = 4.6686e-10*Input->Lines[i].Geff\
          *Input->Lines[i].Lambda0*Input->Lines[i].Lambda0;
    }

    Input->fpixel[0] = 1;
    Input->fpixel[1] = 1;
    Input->fpixel[2] = 1;
    Input->fpixel[3] = 1;  

    Input->npar = 0;
    for(i=1;i<10;i++){
      if(Input->inv[i]) Input->npar++;
    }

    Input->step[1] = 1000;
    Input->step[2] = 2.0;
    Input->step[3] = 2.0;
    Input->step[4] = 2.0;
    Input->step[5] = 20.;
    Input->step[6] = 0.3;
    Input->step[7] = 30;
    Input->step[8] = 1e5;
    Input->step[9] = 0.3;

    fclose(Fa);

    if(Mpi->nprocs==1){
      LM->Hessian = (double **)MATRIX(1, Input->npar, 1, \
          Input->npar, enum_dbl, false);
      LM->Hessian_new = (double **)MATRIX(1, Input->npar, 1, \
          Input->npar, enum_dbl, false);
      LM->Jacfvec = (double *)VECTOR(1, Input->npar, enum_dbl, false);
      LM->Regul_H = (double **)MATRIX(1, Input->npar, 1, \
          Input->npar, enum_dbl, true);
      LM->Regul_J = (double *)VECTOR(1, Input->npar, enum_dbl, true);
      LM->Sol = (double *)VECTOR(1, Input->npar, enum_dbl, false);
      LM->indx = (int *)VECTOR(1, Input->npar, enum_int, false);
      Input->W = (double *)VECTOR(1, Input->npar, enum_dbl, false);
      Input->V = (double **)MATRIX(1, Input->npar, 1, Input->npar, \
          enum_dbl, false);
      Input->Fadd = (STRUCT_FADDEEVA *)malloc(sizeof(STRUCT_FADDEEVA));
      Input->Fadd->nfigure = 6;
      Faddeeva_init(Input->Fadd);
    }else{
      if(Mpi->rank>0){
        LM->Hessian = (double **)MATRIX(1, Input->npar, 1, \
            Input->npar, enum_dbl, false);
        LM->Hessian_new = (double **)MATRIX(1, Input->npar, 1, \
            Input->npar, enum_dbl, false);
        LM->Jacfvec = (double *)VECTOR(1, Input->npar, enum_dbl, false);
        LM->Regul_H = (double **)MATRIX(1, Input->npar, 1, \
            Input->npar, enum_dbl, true);
        LM->Regul_J = (double *)VECTOR(1, Input->npar, enum_dbl, true);
        LM->Sol = (double *)VECTOR(1, Input->npar, enum_dbl, false);
        LM->indx = (int *)VECTOR(1, Input->npar, enum_int, false);
        Input->W = (double *)VECTOR(1, Input->npar, enum_dbl, false);
        Input->V = (double **)MATRIX(1, Input->npar, 1, Input->npar, \
            enum_dbl, false);
        Input->Fadd = (STRUCT_FADDEEVA *)malloc(sizeof(STRUCT_FADDEEVA));
        Input->Fadd->nfigure = 6;
        Faddeeva_init(Input->Fadd);
      }else{
        free(LM);
      }
    }

    return 0;

}

/*--------------------------------------------------------------------------------*/

extern int Free_Ram(STRUCT_LM *LM, STRUCT_STK *Stk, STRUCT_INPUT *Input, 
    STRUCT_PAR *Par, STRUCT_MPI *Mpi){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        Free the memory.
      Record of revisions:
        28 Jun. 2024 (Hao Li)  
      Input parameters:
        LM, structure with the hessian matrix.
        Stk, structure with the Stokes profiles
        Input, structure with the input information.
        Par, structure with the parameter information.
        Mpi, structure with Mpi information.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/
    
    if(Mpi->nprocs==1){
      FREE_VECTOR(LM->Sol, 1, enum_dbl);
      FREE_VECTOR(LM->indx, 1, enum_int);
      FREE_VECTOR(LM->Jacfvec, 1, enum_dbl);
      FREE_MATRIX((void *)LM->Hessian, 1, 1, enum_dbl);
      FREE_MATRIX((void *)LM->Hessian_new, 1, 1, enum_dbl);
      FREE_VECTOR(LM->Regul_J, 1, enum_dbl);
      FREE_MATRIX((void *)LM->Regul_H, 1, 1, enum_dbl);
      free(LM);

      FREE_VECTOR(Input->W, 1, enum_dbl);
      FREE_MATRIX((void *)Input->V, 1, 1, enum_dbl);

      free(Input->Fadd->Expa2n2);
      free(Input->Fadd);

      FREE_TENSOR_DBL(Stk->Jacobian, 0, 0, 0);
      free(Stk->Lambda);
      free(Stk->noise);

      free(Input->data_int);
      free(Input->data_flt);
    }else{
      if(Mpi->rank>0){
        FREE_VECTOR(LM->Sol, 1, enum_dbl);
        FREE_VECTOR(LM->indx, 1, enum_int);
        FREE_VECTOR(LM->Jacfvec, 1, enum_dbl);
        FREE_MATRIX((void *)LM->Hessian, 1, 1, enum_dbl);
        FREE_MATRIX((void *)LM->Hessian_new, 1, 1, enum_dbl);
        FREE_VECTOR(LM->Regul_J, 1, enum_dbl);
        FREE_MATRIX((void *)LM->Regul_H, 1, 1, enum_dbl);
        free(LM);

        FREE_VECTOR(Input->W, 1, enum_dbl);
        FREE_MATRIX((void *)Input->V, 1, 1, enum_dbl);

        free(Input->Fadd->Expa2n2);
        free(Input->Fadd);

        FREE_TENSOR_DBL(Stk->Jacobian, 0, 0, 0);
        free(Stk->Lambda);
        free(Stk->noise);

      }else{
        free(Input->data_int);
        free(Input->data_flt);
      }
    }
    
    FREE_MATRIX((void *)Stk->prof, 0, 0, enum_dbl);
    FREE_MATRIX((void *)Stk->fit, 0, 0, enum_dbl);
    FREE_MATRIX((void *)Stk->syn, 0, 0, enum_dbl);
    free(Stk);

    FREE_VECTOR(Par->Par, 1, enum_dbl);
    FREE_VECTOR(Par->Par_Best, 1, enum_dbl);
    FREE_VECTOR(Par->Par_tmp, 1, enum_dbl);
    
    free(Input->Lines);
    free(Input);
    free(Mpi);

    return 0;

}

/*--------------------------------------------------------------------------------*/
