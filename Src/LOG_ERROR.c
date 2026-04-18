
#include "LOG_ERROR.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:

        02 Mar. 2026  (Hao Li)
          --- Updates: 
              redesign logging and error handling module. 

        24 Apr. 2024
          --- Initial commit (Hao Li)
     
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

static FILE *log_file = NULL;
static char buffer[MAX_BUFRER_SIZE];
char MeSS[MAX_MESSAGE_LENGTH];

/*--------------------------------------------------------------------------------*/

void ABORTED(void){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        terminate the job.
      Record of revisions:
        02 Mar. 2026.
      Input parameters:
        .
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    #ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
    #else
      exit(EXIT_FAILURE);                      
    #endif
}

/*--------------------------------------------------------------------------------*/

void LOG_INIT(const char *filename){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        Caching log file handles.
      Record of revisions:
        02 Mar. 2026.
      Input parameters:
        filename, the log file.
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(!filename) return;

    if(filename){
      log_file = fopen(filename, "a");
      if(!log_file){ 
        log_file = fopen(filename, "w");
        if(!log_file){ 
          LOG_ERROR(ERR_LVL_ERROR, "LOG_INIT", "Failed to open log file");
        }
      }
    }

    return;
}

/*--------------------------------------------------------------------------------*/

void LOG_FINALIZE(void){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        close the log file.
      Record of revisions:
        02 Mar. 2026.
      Input parameters:
        .
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(log_file){
      fclose(log_file);
      log_file = NULL;
    }

    return;
}

/*--------------------------------------------------------------------------------*/

void LOG_WRITE(const char *msg, bool to_screen, bool verbose_flag){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        close the log file.
      Record of revisions:
        02 Mar. 2026.
      Input parameters:
        msg, the log message.
        to_screen, write message to screen or to the log file.
        verbose_flag, verbose mode or not.
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if (verbose_flag) {
      if(to_screen) fprintf(stderr, "%s\n", msg);
      if(log_file) fprintf(log_file, "%s\n", msg);
    }

    return;
}

/*--------------------------------------------------------------------------------*/

void LOG_ERROR(ERR_LVL lvl, const char *rname, const char *msg){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        write error log.
      Record of revisions:
        02 Mar. 2026.
      Input parameters:
        lvl,  error level (enum_warning, enum_error)
        rname, subroutine name. 
        msg, the log message.
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    snprintf(buffer, sizeof(buffer), "%s in routine %s: %s",
             lvl==ERR_LVL_ERROR?"-ERROR":"-WARNING",
             rname, msg);

    LOG_WRITE(buffer, true, true);

    if(lvl==ERR_LVL_ERROR) ABORTED();

    return;
}

/*--------------------------------------------------------------------------------*/

void LOG_MODEL(const double *Par, bool verbose_flag){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        write model parameters.
      Record of revisions:
        02 Mar. 2026.
      Input parameters:
        Par, the model parameter.
        verbose_flag, verbose mode or not.
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    #define Par_Rad2Deg   57.29577951308232286464

    if(!verbose_flag) return;

    LOG_WRITE("       -- model parameters -- ", false, true);
    snprintf(buffer, sizeof(buffer), "  - field strength: %.2f (G)", Par[0]); 
    LOG_WRITE(buffer, false, true);
    snprintf(buffer, sizeof(buffer), 
        "  - inclination: %.2f (deg)", Par[1]*Par_Rad2Deg); 
    LOG_WRITE(buffer, false, true);
    snprintf(buffer, sizeof(buffer), 
        "  - azimuth: %.2f (deg)", Par[2]*Par_Rad2Deg); 
    LOG_WRITE(buffer, false, true);
    snprintf(buffer, sizeof(buffer), "  - velocity: %.3f (km/s)", Par[3]); 
    LOG_WRITE(buffer, false, true);
    snprintf(buffer, sizeof(buffer), "  - Doppler width: %.3f (mA)", Par[4]); 
    LOG_WRITE(buffer, false, true);
    snprintf(buffer, sizeof(buffer), "  - damping width: %.3f (Dopp)", Par[5]); 
    LOG_WRITE(buffer, false, true);
    snprintf(buffer, sizeof(buffer), "  - Eta: %.2f", Par[6]); 
    LOG_WRITE(buffer, false, true);
    snprintf(buffer, sizeof(buffer), "  - continuum: %e", Par[7]); 
    LOG_WRITE(buffer, false, true);
    snprintf(buffer, sizeof(buffer), "  - factor beta: %e", Par[8]); 
    LOG_WRITE(buffer, false, true);

    return;
}

/*--------------------------------------------------------------------------------*/

bool FILE_EXIST(const char *Filename){

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
