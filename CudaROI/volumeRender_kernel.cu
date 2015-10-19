/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

// Simple 3D volume renderer

#ifndef _VOLUMERENDER_KERNEL_CU_
#define _VOLUMERENDER_KERNEL_CU_

#include "reductionMax.hh"
#include "volumeRender_kernel.cuh"

typedef unsigned int  uint;
typedef unsigned char uchar;

typedef struct {
    float4 m[3];
} float3x4;

typedef unsigned short VolumeType;
//typedef float VolumeType;

__constant__ float3x4 c_invViewMatrix;  // inverse view matrix
texture<VolumeType, 3, cudaReadModeElementType> tex;         // 3D texture
texture<VolumeType, 3, cudaReadModeElementType> tex_cluster;         // 3D texture

struct Ray {
	float3 o;	// origin
	float3 d;	// direction
};

// intersect ray with a box
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm

__device__
int intersectBox(Ray r, float3 boxmin, float3 boxmax, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float3 invR = make_float3(1.0f) / r.d;
    float3 tbot = invR * (boxmin - r.o);
    float3 ttop = invR * (boxmax - r.o);

    // re-order intersections to find smallest and largest on each axis
    float3 tmin = fminf(ttop, tbot);
    float3 tmax = fmaxf(ttop, tbot);

    // find the largest tmin and the smallest tmax
    float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), fmaxf(tmin.x, tmin.z));
    float smallest_tmax = fminf(fminf(tmax.x, tmax.y), fminf(tmax.x, tmax.z));

	*tnear = largest_tmin;
	*tfar = smallest_tmax;

	return smallest_tmax > largest_tmin;
}

// transform vector by matrix (no translation)
__device__
float3 mul(const float3x4 &M, const float3 &v)
{
    float3 r;
    r.x = dot(v, make_float3(M.m[0]));
    r.y = dot(v, make_float3(M.m[1]));
    r.z = dot(v, make_float3(M.m[2]));
    return r;
}

// transform vector by matrix with translation
__device__
float4 mul(const float3x4 &M, const float4 &v)
{
    float4 r;
    r.x = dot(v, M.m[0]);
    r.y = dot(v, M.m[1]);
    r.z = dot(v, M.m[2]);
    r.w = 1.0f;
    return r;
}

__device__
float4 color_interpolate_cluster(float sample){
	
	// ACCENT
	
	// if(sample <= 1)
	// 	return make_float4((float)0.99215,(float)0.75294, (float)0.52549, 1.0);
	// else if(sample <= 2)
	// 	return make_float4( (float)0.498, (float)0.7882, (float)0.498, 0.25);
	// else if(sample <= 3)
	// 	return make_float4((float)0.74509,(float)0.68235, (float)0.83137, 1.0);
	// else if(sample <= 4)
	// 	return make_float4(1.0,1.0,1.0,1.0);
	
	// Dark2
	
	if(sample <= 1)
		return make_float4( 0.8509803921569,0.3725490196078,0.007843137254902, 1.0);
	else if(sample <= 2)
		return make_float4( 0.1058823529412, 0.6196078431373, 0.4666666666667, 0.25);
	else if(sample <= 3)
		return make_float4( 0.4588235294118,0.4392156862745,0.7019607843137, 1.0);
	else if(sample <= 4)
		return make_float4(1.0,1.0,1.0,1.0);
	
	return make_float4(0.0,0.0,0.0,0.0);	
}
__device__ 
float4 color_interpolate_large(float sample, float4 one, float4 two, float4 three,
				float4 four, float4 five, float4 six){
	
	float4 retcolor = make_float4(0);
	float percent = 0.0f; 
		
	if(sample <= 0.2f){
	
		percent = (0.2f - sample) / 0.2f;
		retcolor = (percent)*one + (1.0f-percent) * two;
		
	}else if(sample > 0.2f && sample <= 0.3f){
		
		percent = (0.3f - sample)  / 0.1f;
		retcolor = (percent)*two + (1.0f-percent) * three;
		
	}else if(sample > 0.3f && sample <= 0.4f){
		
		percent = (0.4f - sample) / 0.1f;
		retcolor = (percent)*three + (1.0f-percent) * four;
		
	}else if(sample > 0.4f && sample <= 0.5f){
		
		percent = (0.5f - sample) / 0.1f;
		retcolor = (percent)*four + (1.0f-percent) * five;
		
	}else{
		
		percent = (1.0 - sample) / 0.5f;
		retcolor = (percent)*five + (1.0f-percent) * six;
	}
	
	return retcolor;	
}
__device__ 
float4 color_interpolate(float sample, float4 one, float4 two, float4 three,
				float4 four, float4 five, float4 six){
	
	float4 retcolor = make_float4(0);
	float percent = 0.0f; 
		
	if(sample <= 25500.0f){
	
		percent = (25500.0f - sample) / 25500.0f;
		retcolor = (percent)*one + (1.0f-percent) * two;
		
	}else if(sample > 25500.0f && sample <= 26500.0f){
		
		percent = (26500.0f - sample)  / 1000.0f;
		retcolor = (percent)*two + (1.0f-percent) * three;
		
	}else if(sample > 26500.0f && sample <= 27500.0f){
		
		percent = (27500.0f - sample) / 1000.0f;
		retcolor = (percent)*three + (1.0f-percent) * four;
		
	}else if(sample > 27500.0f && sample <= 28500.0f){
		
		percent = (28500.0f - sample) / 1000.0f;
		retcolor = (percent)*four + (1.0f-percent) * five;
		
	}else{
		
		percent = (65535.0f - sample) / 65535.0f;
		retcolor = (percent)*five + (1.0f-percent) * six;
	}
	
	return retcolor;	
}

__device__ uint rgbaFloatToInt(float4 rgba, float global_max, float red, float green, float blue)
{
	rgba.x = rgba.x / (global_max+2);
    rgba.y = rgba.y / (global_max+2);
    rgba.z = rgba.z / (global_max+2);
	rgba.w = 0.5;
	
    return (uint(rgba.w*255)<<24) | (uint(rgba.z*255)<<16) | (uint(rgba.y*255)<<8) | uint(rgba.x*255);
}

__global__ void
d_render(float4 *d_iColors, ushort *data,
 						float *d_iRed, float *d_iGreen, float *d_iBlue, uint imageW, uint imageH,
     					float density, float brightness, float4 one, float4 two, float4 three, 
						float4 four, float4 five, float4 six, int type)
{
    const int maxSteps = 500;
    const float tstep = 0.01f;
    const float3 boxMin = make_float3(-1.0f, -1.0f, -1.0f);
    const float3 boxMax = make_float3(1.0f, 1.0f, 1.0f);

    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;
    if ((x >= imageW) || (y >= imageH)) return;

    float u = (x / (float) imageW)*2.0f-1.0f;
    float v = (y / (float) imageH)*2.0f-1.0f;

    // calculate eye ray in world space
    Ray eyeRay;
    eyeRay.o = make_float3(mul(c_invViewMatrix, make_float4(0.0f, 0.0f, 0.0f, 1.0f)));
    eyeRay.d = normalize(make_float3(u, v, -2.0f));
    eyeRay.d = mul(c_invViewMatrix, eyeRay.d);

    // find intersection with box
	float tnear, tfar;
	int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);
    if (!hit) return;
	if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from front to back, accumulating color
    float4 sum = make_float4(0.0f);
    float t = tnear;
    float3 pos = eyeRay.o + eyeRay.d*tnear;
    float3 step = eyeRay.d*tstep;
    float sample = 0;

    for(int i=0; i<maxSteps; i++) {
        
		// read from 3D texture
       // remap position to [0, 1] coordinates
  		if(type == 0)
			sample = tex3D(tex, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);
		else
			sample = tex3D(tex_cluster, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);
			float4 col = make_float4(0.0f);

        // lookup in transfer function texture
		if(type == 0)
			col = color_interpolate(sample,one,two,three,four,five,six);
		else
			col = color_interpolate_cluster(sample);
        // pre-multiply alpha
		col.x *= col.w;
		col.y *= col.w;
		col.z *= col.w;
        // "over" operator for front-to-back blending
		sum = sum + col;//*(1.0f - sum.w);

        t += tstep;
        if (t > tfar) break;
		
        pos += step;
    }

    sum *= brightness;

    d_iColors[y*imageW + x] = sum;
 
    d_iRed[y*imageW + x] = sum.x;
    d_iGreen[y*imageW + x] = sum.y;
    d_iBlue[y*imageW + x] = sum.z;
}

__global__
void create_image(uint *output, float4 *d_iColors, float global_max, float red, float green, float blue, uint imageW, uint imageH){

    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;
    if ((x >= imageW) || (y >= imageH)) return;

    output[y*imageH+x] = rgbaFloatToInt(d_iColors[y*imageW+x], global_max, red, green, blue);
} 

void setup_cluster(void *cluster, cudaExtent volumeSize, uint image_size, cudaArray *d_volumeArray_cluster){
	
	// Cluster setup

		// create 3D array
		cudaChannelFormatDesc channelDesc_cluster = cudaCreateChannelDesc<VolumeType>();
		cutilSafeCall( cudaMalloc3DArray(&d_volumeArray_cluster, &channelDesc_cluster, volumeSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr(cluster, volumeSize.width*sizeof(VolumeType), volumeSize.width, volumeSize.height);
		copyParams.dstArray = d_volumeArray_cluster;
		copyParams.extent   = volumeSize;
		copyParams.kind     = cudaMemcpyHostToDevice;
		cutilSafeCall( cudaMemcpy3D(&copyParams) );  

		// set texture parameters
		tex_cluster.normalized = true;                      // access with normalized texture coordinates
		tex_cluster.filterMode = cudaFilterModePoint;      // linear interpolation
		tex_cluster.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
		tex_cluster.addressMode[1] = cudaAddressModeClamp;

		// bind array to 3D texture
		cutilSafeCall(cudaBindTextureToArray(tex_cluster, d_volumeArray_cluster, channelDesc_cluster));
	
}

void setup_volume(void *h_volume, cudaExtent volumeSize, uint image_size, cudaArray *d_volumeArray){
	
	// create 3D array
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<VolumeType>();
    cutilSafeCall( cudaMalloc3DArray(&d_volumeArray, &channelDesc, volumeSize) );

    // copy data to 3D array
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr   = make_cudaPitchedPtr(h_volume, volumeSize.width*sizeof(VolumeType), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_volumeArray;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    cutilSafeCall( cudaMemcpy3D(&copyParams) );  

    // set texture parameters
    tex.normalized = true;                      // access with normalized texture coordinates
    tex.filterMode = cudaFilterModePoint;      // linear interpolation
    tex.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    tex.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    cutilSafeCall(cudaBindTextureToArray(tex, d_volumeArray, channelDesc));
}

void render_kernel(dim3 gridSize, dim3 blockSize, uint *d_output, uint *d_cluster, float* d_iRed, float* d_oRed, 
						float* d_iGreen, float* d_oGreen, float* d_iBlue, float* d_oBlue, float4* d_iColors, unsigned short* data, 
						unsigned short *cluster_data, uint imageW, uint imageH, float density, float brightness, 
						float4 one, float4 two, float4 three, float4 four, float4 five, float4 six,
						void *h_volume, void *cluster, cudaExtent volumeSize, cudaArray *d_volumeArray, cudaArray *d_volumeArray_cluster, int *set)
{	
		
	int size = imageH * imageW;
	
	if(set[0] == 0){
		setup_volume(h_volume, volumeSize, size, d_volumeArray);
		set[0] = 1;
	}
	if(set[1] == 0){
		setup_cluster(cluster, volumeSize, size, d_volumeArray_cluster);
		set[1] = 1;
	}
	/* clear colors buffers */

	cutilSafeCall(cudaMemset(d_iColors, 0, imageH*imageW*sizeof(float4)));	
	
	cutilSafeCall(cudaMemset(d_iRed, 0, imageH*imageW*sizeof(float)));	
	cutilSafeCall(cudaMemset(d_oRed, 0, imageH*imageW*sizeof(float)));	
	cutilSafeCall(cudaMemset(d_iGreen, 0, imageH*imageW*sizeof(float)));	
	cutilSafeCall(cudaMemset(d_oGreen, 0, imageH*imageW*sizeof(float)));
	cutilSafeCall(cudaMemset(d_iBlue, 0, imageH*imageW*sizeof(float)));	
	cutilSafeCall(cudaMemset(d_oBlue, 0, imageH*imageW*sizeof(float)));

	d_render<<<gridSize, blockSize>>>(d_iColors, data, d_iRed, d_iGreen, d_iBlue, imageW, imageH, density, brightness, 
						one, two, three, four, five, six, 0);

	float max_red = reduce_max(d_oRed, d_iRed, size);
	float max_green = reduce_max(d_oGreen, d_iGreen, size);
	float max_blue = reduce_max(d_oBlue, d_iBlue, size);
	
	float global_max = fmax(max_red, max_green);
	global_max = fmax(global_max, max_blue);
	
	create_image<<<gridSize, blockSize>>>(d_output, d_iColors, global_max, max_red, max_green, max_blue, imageW, imageH);
	
	// render image
	// 
	d_render<<<gridSize, blockSize>>>(d_iColors, cluster_data, d_iRed, d_iGreen, d_iBlue, imageW, imageH, density, brightness, 
					one, two, three, four, five, six, 1);

	max_red = reduce_max(d_oRed, d_iRed, size);
	max_green = reduce_max(d_oGreen, d_iGreen, size);
	max_blue = reduce_max(d_oBlue, d_iBlue, size);

	global_max = fmax(max_red, max_green);
	global_max = fmax(global_max, max_blue);

	create_image<<<gridSize, blockSize>>>(d_cluster, d_iColors, global_max, max_red, max_green, max_blue, imageW, imageH);	  
}
void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix)
{
    cutilSafeCall( cudaMemcpyToSymbol(c_invViewMatrix, invViewMatrix, sizeofMatrix) );
}

#endif // #ifndef _VOLUMERENDER_KERNEL_CU_
