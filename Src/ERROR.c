
#include "ERROR.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:

        24 Apr. 2024
          --- Initial commit (Hao Li)
     
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

char MeSS[MAX_MESSAGE_LENGTH];

/*--------------------------------------------------------------------------------*/

extern int Error(enum error_level error_lv, const char *rname,\
        const char *Emes, char *Filename){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        error handler.
      Record of revisions:
        16 Jan. 2024.
      Input parameters:
        error_lv, level of the error (enum_warning, enum_error)
        rname, name of the routine.
        Emes, the message.
      Return:
        .
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    FILE *Fa;

    switch(error_lv){ 

      case enum_error:
        if(Filename){
          Fa = fopen(Filename, "w+");
          if(!Fa) Fa = fopen(Filename, "w");
          fprintf(Fa, "\n-ERROR in routine %s\n\n %s\n", \
              rname, Emes);
          fclose(Fa);
        }else{
          fprintf(stderr, "\n-ERROR in routine %s\n\n %s\n", \
            rname, Emes);
        }
        ABORTED();
        break;
      case enum_warning:
        if(Filename){
          Fa = fopen(Filename, "w+");
          if(!Fa) Fa = fopen(Filename, "w");
          fprintf(stderr, "\n-WARNING in routine %s\n %s\n", rname, \
              Emes);
          fclose(Fa);
        }else{
          fprintf(stderr, "\n-WARNING in routine %s\n %s\n", rname, \
              Emes);
        }

        return 1;
        
      default:
        return 0;
    }
    return 0;

}

/*--------------------------------------------------------------------------------*/

extern int VerboseM(const char *Emes, char *Filename, bool screen){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        verbose on screen or in a file.
      Record of revisions:
        16 Jan. 2024.
      Input parameters:
        Emes, the message.
        Filename, file name.
        screen, ture on the screen, false in a file.
      Return:
        .
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(screen){

      fprintf(stderr, " %s\n ", Emes);
    }else{

      FILE *Fa = fopen(Filename, "a");

      if(!Fa) Fa = fopen(Filename, "w");

      fprintf(Fa, " %s\n ", Emes);

      fclose(Fa);
    }
  
    return 0;

}

/*--------------------------------------------------------------------------------*/

extern int Verbose(const char *Emes, char *Filename, bool cut){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        verbose in a file or not.
      Record of revisions:
        16 Jan. 2024.
      Input parameters:
        Emes, the message.
        Filename, file name.
        cut, ture verbose, false not.
      Return:
        .
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(!cut){

      FILE *Fa = fopen(Filename, "a");

      if(!Fa) Fa = fopen(Filename, "w");

      fprintf(Fa, " %s\n ", Emes);

      fclose(Fa);
    }
  
    return 0;

}

/*--------------------------------------------------------------------------------*/

extern int Verbose_model(double *Par, char *Filename, bool cut){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        verbose the model.
      Record of revisions:
        16 Jan. 2024.
      Input parameters:
        Par, the model parameter.
        Filename, file name.
        cut, ture verbose, false not.
      Return:
        .
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(!cut){

      FILE *Fa = fopen(Filename, "a");

      if(!Fa) Fa = fopen(Filename, "w");

      fprintf(Fa, "     -- the model parameter -- \n\n");
      fprintf(Fa, "  - field strength : %.2f (G) \n", Par[1]);
      fprintf(Fa, "  - inclination : %.2f (deg) \n", Par[2]/Par_Pi*180.);
      fprintf(Fa, "  - azimuth : %.2f (deg) \n", Par[3]/Par_Pi*180.);
      fprintf(Fa, "  - velocity : %.3f (km/s) \n", Par[4]);
      fprintf(Fa, "  - Doppler width : %.3f (mA) \n", Par[5]);
      fprintf(Fa, "  - damping width : %.3f (Dopp) \n", Par[6]);
      fprintf(Fa, "  - Eta : %2f \n", Par[7]);
      fprintf(Fa, "  - continuum : %e \n", Par[8]);
      fprintf(Fa, "  - factor beta : %e \n\n", Par[9]);

      fclose(Fa);
    }
  
    return 0;

}

/*--------------------------------------------------------------------------------*/

extern bool FILE_EXIST(char *Filename){

/*--------------------------------------------------------------------------------*/
      
    /*######################################################################
      Purpose:
        check if a file exists.
      Record of revisions:
        23 Apr. 2024.
      Input parameters:
        Filename, name of the file
      Return:
        true or false
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    FILE *Fa = fopen(Filename, "r");
    
    if(Fa){

      fclose(Fa);

      return true;
    }

    return false;
}

/*--------------------------------------------------------------------------------*/

