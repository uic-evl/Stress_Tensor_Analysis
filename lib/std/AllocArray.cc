// AllocArray.c

// Routines to allocate arrays in memory dynamically.

// This package also includes numerical formatting routine(s)

#include "AllocArray.hh"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

void** AllocArray2D(int nRows,
                    int nCols,
                    int nSize)
{
    void**  vppArray;
    char*   cpData;
    size_t  iRowTableSize;
    size_t  iDataArraySize;

    size_t  i;

    // calculate the size of the pieces
    iRowTableSize   = (size_t)nRows * sizeof(void*);
    iDataArraySize  = (size_t)nRows * (size_t)nCols * (size_t)nSize + 64;   // allow for cache alignment

    // allocate the memory
    vppArray = (void**)calloc(iRowTableSize + iDataArraySize, 1);

    // confirm allocation
    if (vppArray == NULL)
    {
        fprintf(stderr,
#ifdef _X64_COMPILER
                "\n***** Unable to allocate memory block of %I64u bytes.\n",
#else
                "\n***** Unable to allocate memory block of %I32lu bytes.\n",
#endif
                iRowTableSize + iDataArraySize);
        QUIT;
    }

    // set up the row table
    cpData  = (char*)vppArray;
    cpData += iRowTableSize + 63;
#ifdef _X64_COMPILER
    cpData  = (char*)(((__int64)cpData) & 0xffffffffffffffc0);  // round back to 64-byte boundary
#else
    cpData  = (char*)(((uintptr_t)cpData) & 0xffffffc0);     // round back to 64-byte boundary
#endif

    // calculate the pointers
    for (i=0; i<(size_t)nRows; i++)
    {
        vppArray[i] = &cpData[i * nCols * nSize];
    }

    return vppArray;
}


void*** AllocArray3D(int nPlanes,
                     int nRows,
                     int nCols,
                     int nSize)
{
    void***         vpppArray;
    void**          vppRowTable;
    char*           cpData;
    size_t          iPlaneTableSize;
    size_t          iRowTableSize;
    size_t          iDataArraySize;
    size_t          i,j;
#ifndef _X64_COMPILER
    double          dMaxAllocSize;
#endif

    // calculate the size of the pieces
    iPlaneTableSize = (size_t)nPlanes   * 
                      sizeof(void**);
    iRowTableSize   = (size_t)nPlanes   * 
                      (size_t)nRows     * 
                      sizeof(void*);
    iDataArraySize  = (size_t)nPlanes   * 
                      (size_t)nRows     *    
                      (size_t)nCols     * 
                      (size_t)nSize     + 
                      64;               // allow for cache alignment

#ifndef _X64_COMPILER
    // check for the maximum single-allocation block size in 32-bit machine
    dMaxAllocSize   = (double)nPlanes *
                      (double)nRows   *
                      (double)nCols   *
                      (double)nSize   +
                      64              +
                      iPlaneTableSize +
                      iRowTableSize   ;
    if (dMaxAllocSize > (pow(2.0, 32.0) - 1.0))
    {
        // This array is too large for a single allocation in 32-bit Windows.  
        // The number of bytes will overflow a 32-bit unsigned integer!
        fprintf(stderr, "\n***** 3D ARRAY SIZE EXCEEDS 2^32 BYTES. *****\n");
        QUIT;
    }
#endif

    // allocate the memory 
    vpppArray = (void***)calloc(iPlaneTableSize + iRowTableSize + iDataArraySize, 1);

    // confirm allocation
    if (vpppArray == NULL)
    {
#ifdef _X64_COMPILER
        fprintf(stderr,
                "\n***** Unable to allocate memory block of %I64u bytes.\n",
                iPlaneTableSize + iRowTableSize + iDataArraySize);
        QUIT;
#else   // #ifdef _X64_COMPILER
        // allocation failed!  Let's try it again in large-array mode.
        goto AlternativeAllocationMethod;
#endif
    }

    // populate the plane table
    vppRowTable = (void**) &vpppArray[nPlanes];             // start just past the plane table
    cpData      = (char*)  &vppRowTable[nPlanes * nRows];   // start just past the row table

	// align the start of the data table with the cache boundaries (round-up to a multiple of 64)
    cpData += 63;
#ifdef _X64_COMPILER
    cpData  = (char*)(((__int64)cpData) & 0xffffffffffffffc0);  // round back to 64-byte boundary
	printf("64 bit son!\n");
#else
    cpData  = (char*)(((uintptr_t)cpData) & 0xffffffc0);     // round back to 64-byte boundary
#endif

    // now calculate the pointers
    for (j=0; j<(size_t)nPlanes; j++)
    {
        vpppArray[j] = &vppRowTable[j * nRows];
        // now populate the row table
        for (i=0; i<(size_t)nRows; i++)
        {
            vpppArray[j][i] = &cpData[(i * nCols + j * nRows * nCols) * nSize];
        }
    }

    return vpppArray;

#ifndef _X64_COMPILER
AlternativeAllocationMethod:
    {
        // This array is so large that we must allocate it plane by plane.
        // The entire 3D array exceeds the capability of 32-bit Windows for a single allocation.
        fprintf(stderr,
                "\nSWITCHING TO ALTERNATIVE ALLOCATION MODE FOR 3D ARRAY OF %.0f BYTES.\n"
                "NOTE THAT THIS ARRAY CAN NEVER BE RELEASED!\n\n",
                dMaxAllocSize);
        // allocate the plane table
        vpppArray = (void***)calloc(nPlanes,sizeof(void**));
        // populate the plane table
        for (j=0; j<(size_t)nPlanes; j++)
        {
            vpppArray[j] = AllocArray2D(nRows, nCols, nSize);
        }
        // we are done!
        return vpppArray;
    }
#endif
}



// convert an integer into a comma-separated ASCII string
char* Comma(int iValue)
{
#define Comma_NumberOfBuffers   20
#define Comma_SizeOfBuffer      20
    static char szBuffer[Comma_NumberOfBuffers][Comma_SizeOfBuffer];
    static int  iBfr = 0;
    char        szWork[Comma_SizeOfBuffer];
    char*       szpOut;
    int         iSrc,iDest;
    // index the buffer selector
    {
        iBfr++;
        if (iBfr >= Comma_NumberOfBuffers) iBfr = 0;
        szpOut = szBuffer[iBfr];
    }
    // convert the radix
    sprintf(szWork,
            "%15d",
            iValue);
    // format the output string
    for (iSrc=iDest=0; szWork[iSrc]; iSrc++)
    {
        switch (szWork[iSrc])
        {
        case ' ':
            break;
        case '-':
            szpOut[iDest++] = '-';
            break;
        default:
            szpOut[iDest++] = szWork[iSrc];
            if ((iSrc%3) == 2)
                szpOut[iDest++] = ',';
            break;
        }
    }
    szpOut[iDest-1] = 0;
    // return a pointer to the buffer
    return szpOut;
}
