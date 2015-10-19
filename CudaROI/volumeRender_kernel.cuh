#ifndef _VOLUMERENDER_KERNEL_H_
#define _VOLUMERENDER_KERNEL_

#include <cutil_inline.h>    // includes cuda.h and cuda_runtime_api.h
#include <cutil_math.h>
#include <stdio.h>

void render_kernel(dim3 gridSize, dim3 blockSize, uint *d_output, uint *d_cluster,float* d_iRed, float* d_oRed, 
						float* d_iGreen, float* d_oGreen, float* d_iBlue, float* d_oBlue, float4* d_iColors, unsigned short* data, 
						unsigned short* cluster_data, uint imageW, uint imageH, float density, float brightness, 
						float4 one, float4 two, float4 three, float4 four, float4 five, float4 six,
						void *h_volume, void *cluster, cudaExtent volumeSize, cudaArray *d_volumeArray, cudaArray *d_volumeArray_cluster, int *set);

void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix);
#endif // VOLUMERENDER_KERNEL