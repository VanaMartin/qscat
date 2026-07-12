#include<stdio.h>
#include<iostream>
#include<cstring>

#include "sec_stream.h"

void fopen_s(FILE** pFile, const char * filename, const char * mode)
{
    *pFile = fopen(filename, mode);
}

void strcpy_s(char * destination, int num, const char * source)
{
    strcpy(destination, source);
}
