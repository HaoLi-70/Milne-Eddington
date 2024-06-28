
#ifndef STR_h
#define STR_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>

/*--------------------------------------------------------------------------------*/

#define Key_Length 50
#define Max_Line_Length 300

extern int Check_Char(char *str, char c);

extern int Indx_Char(char *str, char c, int order);

extern int Indx_Space(char *str, int order);

extern void Trim1(char *str);

extern void Trim2(char *str);

extern void Trim(char *str, int indx);

extern void String_Copy(char *str1, char *str2, long len, bool trim_flag);

extern int String_Split(char *str1, char *str2);

extern void String_to_Upper(char *str);

extern int String_elements(char *str);

extern int Read_line(char *lines, FILE *fa);

/*--------------------------------------------------------------------------------*/

#endif /* STR_h */
