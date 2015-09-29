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
// 
//  Shaders.cc
//  Project4
//  
//  Created by Timothy Luciani on 2011-03-02.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
// 
 
#define STRINGIFY(A) #A

// vertex shader
const char *vertexShader = STRINGIFY(

void main(void)
{
    // calculate window-space point size

	gl_PointSize = 2.0; // apriore 

    gl_TexCoord[0] = gl_MultiTexCoord0;
    gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);

    gl_FrontColor = gl_Color;

}
);

const char *linePixelShader = STRINGIFY(

	uniform vec4 color;

void main(void)
{		
	const vec3 lightDir1 = vec3(0.577,0.577,0.577);
	const vec3 lightDir2 = vec3(0.577,0.577,-0.577);
	const vec3 lightDir3 = vec3(0.577,-0.577,-0.577);
	const vec3 lightDir4 = vec3(-0.577,-0.577,-0.577);
	const vec3 lightDir5 = vec3(-0.577,-0.577,0.577);
	const vec3 lightDir6 = vec3(-0.577,0.577,0.577);
	const vec3 lightDir7 = vec3(-0.577,-0.577,-0.577);
	
	vec3 N;
   
	// calculate normal from texture coordinates
    N.xy = gl_TexCoord[0].xy*vec2(2.0, -2.0) + vec2(-1.0, 1.0);
    float mag = dot(N.xy, N.xy);
    if (mag > 1.0) discard;   // kill pixels outside circle
  
	N.z = sqrt(1.0-mag);

	float diffuse;
	// calculate lighting
	diffuse = max(0.0, dot(lightDir1, N));
	diffuse += max(0.0, dot(lightDir2, N));
	diffuse += max(0.0, dot(lightDir3, N));
	diffuse += max(0.0, dot(lightDir4, N));
	diffuse += max(0.0, dot(lightDir5, N));
	diffuse += max(0.0, dot(lightDir6, N));
	diffuse += max(0.0, dot(lightDir7, N));
	
	gl_FragColor = color * diffuse;
	
}

);
