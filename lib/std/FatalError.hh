#pragma once

// macro to use FatalError() more easily
#define FATAL_ERROR(x)  FatalError(x,__LINE__,__FILE__)

// declaration of the global FatalError() routine
void FatalError(const char* szpErrorMessage, int iLineNumber, char* szpFileName);
