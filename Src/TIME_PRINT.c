
#include "TIME_PRINT.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        8 Sept. 2021
          --- Initial commit (Hao Li)
     
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

extern int Time_Print(void){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        compute and print the running time.
      Record of revisions:
        8 Sept. 2021 (Hao Li)
      Input parameters:
        .
      Return:
        return the current conts.
    ######################################################################*/
    
/*--------------------------------------------------------------------------------*/

    static int conts = 0;
    clock_t time_tmp;
    static clock_t time_begin = 0;
    
    if(conts == 0){
      time_begin = clock();
      fprintf(stderr, "\n Time Calculating Is Initialized \n");
      conts++;
        
    }else{
      time_tmp = clock();
      long hours = (time_tmp-time_begin)/CLOCKS_PER_SEC/3600;      
      long minutes = ((time_tmp-time_begin)/CLOCKS_PER_SEC-hours*3600)/60;
      double seconds = (time_tmp-time_begin)*1.0 \
          /CLOCKS_PER_SEC-hours*3600-minutes*60;
        
      fprintf(stderr, "\n Time Print Point %d \n",conts);
        
      if(hours > 0){
        fprintf(stderr," Running time= %lu h %lu min %.2lf sec\n", \
            hours,minutes,seconds);
      }else if(minutes > 0){
        fprintf(stderr," Running time= %lu min %.2lf sec\n",minutes,seconds);
      }else{
        fprintf(stderr," Running time= %.2lf sec\n",seconds);
      }
      conts++;
    }
    return conts-1;
}

/*--------------------------------------------------------------------------------*/

