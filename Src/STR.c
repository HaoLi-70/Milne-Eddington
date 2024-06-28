
#include "STR.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        30 Oct. 2022
          --- Initial commit (Hao Li)
     
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

extern int Check_Char(char *str, char c){

/*--------------------------------------------------------------------------------*/
    
    /*######################################################################
      Purpose:
        check how many char c is in a string.
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        str, the input string.
        c, the charaster.
      Return:
        return the position of the charaster.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(str == NULL || *str == '\0') return -1;
    
    char *p = str;
    int indx = 0, num = 0;
    
    while(*p != '\0'){
      if(*p == c) num++;
      p++;
      indx++;
    }
    
    if((int)strlen(str) > indx){
      return num;
    }else{
      return -2;
    }
}

/*--------------------------------------------------------------------------------*/

extern int Indx_Char(char *str, char c, int order){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        get the index of char c in a string.
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        str, the input string.
        c, the charaster.
        order, the order of the charaster.
      Return:
        return the position of the charaster.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(str == NULL || *str == '\0') return -1;

    if(order <= 0) return -2;
    
    char *p = str;
    int indx = 0, num = 0;

    while(*p != '\0' && num!=order ){
      if(*p == c) num++;
      p++;
      indx++;
    }

    if((int)strlen(str) > indx){
      return indx;
    }else{
      return -3;
    }
}

/*--------------------------------------------------------------------------------*/

extern int Indx_Space(char *str, int order){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        get the index of a space in a string.
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        str, the input string.
        order, the order of the space.
      Return:
        return the position of the charaster.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(str == NULL || *str == '\0') return -1;

    if(order <= 0) return -2;
    
    char *p = str;
    int indx = 0, num = 0;

    while(*p != '\0' && num!=order ){
      if(isspace(*p)) num++;
      p++;
      indx++;
    }

    if((int)strlen(str) > indx){
      return indx;
    }else{  
      return -3;
    }
}

/*--------------------------------------------------------------------------------*/

extern void Trim1(char *str){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        remove leading spaces in a string.
      Record of revisions:
        29 Nov. 2019.
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
      Return:
        .
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(str == NULL || *str == '\0') return;
    
    int len = 0;
    char *p = str;
    
    while(*p != '\0' && isspace(*p)){
      p++;
      len++;
    }
    
    memmove(str, p, strlen(str) - len + 1);
    
    return;
}

/*--------------------------------------------------------------------------------*/

extern void Trim2(char *str){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        remove trailing spaces in a string.
      Record of revisions:
        29 Nov. 2019.
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    if(str == NULL || *str == '\0') return;
    
    long len = strlen(str);
    char *p = str + len - 1;

    while(p >= str  && isspace(*p)){
      *p = '\0';
      p--;
    }
    
    return;
}

/*--------------------------------------------------------------------------------*/

extern void Trim(char *str, int indx){

/*--------------------------------------------------------------------------------*/
 
    /*######################################################################
      Purpose:
        remove spaces in a string.
      Record of revisions:
        29 Nov. 2019.
      Input parameters:
        str, the input string.
        indx, indx = 1, for the leading spaces; indx =2, for the
            trailing spaces; indx=3 for both.
      Output parameters:
        str, the output string.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    switch(indx){
      case 1:
        Trim1(str);
        break;
            
      case 2:
        Trim2(str);
        break;
            
      case 3:
        Trim1(str);
        Trim2(str);
        break;
            
      default:
        break;
    }
    
    return;
}

/*--------------------------------------------------------------------------------*/

extern void String_Copy(char *str1, char *str2, long len, bool trim_flag){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        copy a string from str2 to str1.
      Record of revisions:
        29 Nov. 2019 (Hao Li)
      Input parameters:
        str2, the input string.
        len, the length.
        trim_flag, if the flag > 0, the leading and trailing spaces in
            str1 will be removed.
      Input parameters:
        str1, the output string.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    memmove(str1, str2, len);
    str1[len] = '\0';
    
    if(trim_flag) Trim(str1, 3);
    
    return;
}

/*--------------------------------------------------------------------------------*/

extern int String_Split(char *str1, char *str2){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        copy the first element in str2 to str1, the element in str2 is 
            deleted.
      Record of revisions:
        29 Nov. 2019 (Hao Li)
      Input parameters:
        str2, the input string.
      Ouput parameters:
        str1, the coppyed element.
        str2, the left elements.
      Return:
        return -1 if no element left, else return 0.
    ######################################################################*/
    
/*--------------------------------------------------------------------------------*/

    long len_tot;
    int len;
    char *p;
    
    Trim(str2, 3);
    len_tot = strlen(str2);
    len = Indx_Space(str2, 1);

    if(len == -3){
      String_Copy(str1, str2, len_tot, 1);
      str2[0] = '\0';
      return -1;

    }else{
      String_Copy(str1, str2, len-1, 1);

      p = str2+len;
      String_Copy(str2, p, len_tot-len, 1);
      Trim(str2, 3);
      Trim(str1, 3);
    }
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern void String_to_Upper(char *str){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        converts all the charactors in str to their uppercases.
      Record of revisions:
        29 Nov. 2019 (Hao Li)
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
      Return:
        .
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    int i;

    for(i = 0; i < (int)strlen(str); i++) str[i] = toupper(str[i]);

    return;
}

/*--------------------------------------------------------------------------------*/

extern int String_elements(char *str){

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        check how many elements (seperated by space) are the in a string.
      Record of revisions:
        29 Nov. 2019 (Hao Li)
      Input parameters:
        str, the input string.
      Return:
        the number of the elements.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    Trim(str, 3);

    int Num = 1, i;

    if(strlen(str) > 0){
      for(i = 0;  i < (int) strlen(str)-1;  i++){
        if(str[i] == ' ' && str[i+1] != ' ') Num++;
      }
      return Num;
    }else{
      return 0;
    }
}

/*--------------------------------------------------------------------------------*/

extern int Read_line(char *lines, FILE *fa){
    
/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Purpose:
        read a not commented line in a file fa.
      Record of revisions:
        29 Nov. 2019 (Hao Li)
      Input parameters:
        fa, the file.
      Output parameters:
        lines, the output string.
      Return:
        return the length of the line.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

    int len_tot, len;

    do{
      if(fgets(lines, 300, fa) == NULL) return -1;
      Trim(lines, 3);
      len = Indx_Char(lines, '!', 1);
      if(len > 0){
        lines[len-1] = '\0';
        Trim(lines, 3);
      }
      len_tot = (int)(strlen(lines));
    }while((lines[0] == '#') || (lines[0] == '!') || \
        (lines[0] == '*') || (len_tot == 0));

    return len_tot;
}

/*--------------------------------------------------------------------------------*/

