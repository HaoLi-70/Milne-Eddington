
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>

/*--------------------------------------------------------------------------------*/

#define Key_Length 64
#define Max_Line_Length 512

/*--------------------------------------------------------------------------------*/

extern int STR_COUNT_CHAR(const char *str, const char c);

extern int STR_INDEX_CHAR(const char *str, const char c, int order);

extern int STR_INDEX_SPACE(const char *str, int order);

extern int STR_TRIM_LEFT(char *str);

extern int STR_TRIM_RIGHT(char *str);

extern int STR_TRIM(char *str);

extern void STR_COPY(char *dest, size_t destsize, const char *src, \
    size_t srcsize, bool trim_flag);

extern int STR_SPLIT(char *dest, size_t destsize, char *src);

extern void STR_TOUPPER(char *str);

extern int STR_ELEMENTS(char *str);

extern int STR_READ_LINE(char **line, size_t *size, FILE *fa);

/*--------------------------------------------------------------------------------*/


