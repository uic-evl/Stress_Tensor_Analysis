// FatalError.cpp : implementation file
//

#include "FatalError.hh"
#include "assert.h"
#include "stdio.h"

// global routine to elicit a fatal error message box
void FatalError(const char* szpErrorMessage, int iLineNumber, char* szpFileName)
{
    printf("Fatal error: %s\n", szpErrorMessage);
}

