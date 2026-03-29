
#include "STR.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        06 Mar. 2026  (Hao Li)
          --- Updates:  
              Redesign the subroutines to improve safety. 
              STR_READ_LINE now supports reading long very line.  
              by dynamically reallocating the array.

        4 Mar. 2026.
          --- Update: Removed unnecessary variables and redundant 
                      computation, and standardize the names 
                      of subroutines. (Hao Li)
                      Redesigned str_read_line to handle dynamically 
                      allocated buffers. (Hao Li)

        30 Oct. 2022
          --- Initial commit (Hao Li)
     
    ######################################################################*/

/*--------------------------------------------------------------------------------*/


int STR_COUNT_CHAR(const char *str, const char c){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        count char c in a string.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        str, the input string.
        c, the charaster.
      Return:
        return the position of the charaster.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(!str) return -1;
    
    int count = 0;
    const char *p = str;

    while(*p){
      if(*p == c) count++;
      p++;
    }
    
    return count;
}

/*--------------------------------------------------------------------------------*/

int STR_INDEX_CHAR(const char *str, const char c, int order){
    
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        get the index of char c in a string.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        str, the input string.
        c, the charaster.
        order, the order of the charaster.
      Return:
        return the position of the charaster.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(!str || *str=='\0') return -1;
    if(order <= 0) return -2;
    
    int count = 0;
    const char *p = str;
    while(*p){
      if(*p == c && ++count == order)
        return p - str;
      p++;
    }

    return -3;
}

/*--------------------------------------------------------------------------------*/


int STR_INDEX_SPACE(const char *str, int order){
    
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        get the index of a space in a string.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        str, the input string.
        order, the order of the space.
      Return:
        return the position of the charaster.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(str == NULL || !*str) return -1;
    if(order <= 0) return -2;

    int count = 0;
    const char *p = str;
    while(*p){
      if(isspace((unsigned char)*p) && ++count == order)
        return p - str;
      p++;
    }

    return -3;
}

/*--------------------------------------------------------------------------------*/

int STR_TRIM_LEFT(char *str){
    
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        remove leading spaces in a string.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
      Return:
        return the number of removed spaces.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(!str || !*str) return 0;
    
    char *p = str;
    while(*p && isspace((unsigned char)*p)){
      p++;
    }

    int count = p - str; 
    if(count>0){
      memmove(str, p, strlen(p)+1);
    }

    return count;
}

/*--------------------------------------------------------------------------------*/


int STR_TRIM_RIGHT(char *str){
    
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        remove trailing spaces in a string.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
      Return:
        return the number of removed spaces.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(!str || *str=='\0') return 0;
    
   
    char *p_end = str+strlen(str)-1;
    char *p = p_end;
 
    while(p >= str && isspace((unsigned char)*p)){
      p--;
    }
    *(p+1) = '\0';
    
    return (int)(p_end-p);
}

/*--------------------------------------------------------------------------------*/


int STR_TRIM(char *str){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        remove both leading and trailing spaces.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
      Return:
        return the number of removed spaces.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if (!str || !*str) return 0;

    size_t len = strlen(str);
    char *p_start = str;
    char *p_end = str+len-1;
    char *p = p_end;

    while (*p_start && isspace((unsigned char)*p_start)) {
      p_start++;
    }

    while(p >= str && isspace((unsigned char)*p)){
      p--;
    }
    *(p+1) = '\0';
   
    int count_left = p_start-str; 
    int count_right = (int)(p_end-p);

    if (count_left>0) {
      memmove(str, p_start, len-count_left-count_right+1);
    }
    
    return count_left+count_right;
}

/*--------------------------------------------------------------------------------*/


void STR_COPY(char *dest, size_t destsize, const char *src, \
    size_t srcsize, bool trim_flag){
    
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        copy the contents of src to dest.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        destsize, buffer size
        src, the input string.
        srcsize, the src size.
        trim_flag, if the flag > 0, remove the leading and trailing spaces.
      Input parameters:
        dest, the output string.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(!dest || !src || !destsize || !srcsize) return;

    size_t size = destsize>srcsize? srcsize : destsize-1;

    memmove(dest, src, size);

    dest[size] = '\0';
    
    if(trim_flag) STR_TRIM(dest);
    
    return;
}

/*--------------------------------------------------------------------------------*/


int STR_SPLIT(char *dest, size_t destsize, char *src){
    
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        copy the first element in src to dest and remove the copied element 
            from src.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        src, the input string.
        destsize, buffer size.
      Ouput parameters:
        dest, the coppyed element.
        src, the left elements.
      Return:
        return 0 if no element left, else return 1.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(!dest || !src) return -1;

    STR_TRIM(src);
    size_t len_tot = strlen(src);

    int space_idx = STR_INDEX_SPACE(src, 1); 
    if(space_idx < 0){  
      STR_COPY(dest, destsize, src, len_tot, true);
      src[0] = '\0';
      return 0;

    }else{
      STR_COPY(dest, destsize, src, space_idx, true);
      STR_COPY(src, destsize, src + space_idx + 1, len_tot - space_idx - 1, true);
    }

    return 1;
}

/*--------------------------------------------------------------------------------*/


void STR_TOUPPER(char *str){
    
/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        converts all the charactors in str to their uppercases.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
      Return:
        .
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(!str) return; 

    for(char *p=str; *p; p++){
      *p = toupper((unsigned char)*p);  
    }

    return;
}

/*--------------------------------------------------------------------------------*/


int STR_ELEMENTS(char *str){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        check how many elements (seperated by space) are the in a string.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        str, the input string.
      Return:
        the number of the elements.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(!str) return -1; 

    STR_TRIM(str);

    if (!(*str)) return 0;

    int count = 1;
    for(char *p=str;  *(p+1);  p++){
      if(isspace((unsigned char)*p) && !isspace((unsigned char)*(p+1))){ 
        count++;
      }
    }

    return count; 
}

/*--------------------------------------------------------------------------------*/


int STR_READ_LINE(char **line, size_t *size, FILE *fa){

/*--------------------------------------------------------------------------------*/
    /*######################################################################
      Purpose:
        read a valid line from a file fa.
      Record of revisions:
        4 Mar. 2026.
      Input parameters:
        fa, the file.
      Output parameters:
        line, the output string.
        size, bufer size.
      Return:
        return the length of the line.
    ######################################################################*/
/*--------------------------------------------------------------------------------*/

    if(!line || !size || !fa) return -3;

    if(*line == NULL || *size == 0){
      *size = 512;
      *line = malloc(*size);
      if (!*line) return -2;
    }

    while(1){

      size_t len = 0;
      int indexspace, nspaces;

      while(1){
        if(!fgets(*line + len, (int)(*size - len), fa)){
          if(len == 0) return -1; 
          break;
        }

        size_t chunk = strlen(*line+len);
        len += chunk;
        if (len > 0 && (*line)[len-1] == '\n') break;

        size_t new_size = (*size)*2;
        char *tmp = realloc(*line, new_size);
        if (!tmp) return -2;
        *line = tmp;
        *size = new_size;
      }

      if(len == 0) continue;

      nspaces = STR_TRIM_LEFT(*line);
      len -= nspaces;

      char first = (*line)[0];
      if(first=='#' || first=='!' || first=='*' || first=='\0'){
        continue;
      }

      indexspace = STR_INDEX_CHAR(*line, '!', 1);
      if(indexspace >= 0){
        line[indexspace] = '\0';
        if (len>(size_t)indexspace){
          len = indexspace;
        }
      }
      nspaces = STR_TRIM_RIGHT(*line);
      len -= nspaces;

      return (int)len;
    }
}

/*--------------------------------------------------------------------------------*/

