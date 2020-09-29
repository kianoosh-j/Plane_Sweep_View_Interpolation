#ifndef HEADERS
#define HEADERS


#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include "CL/cl.hpp"
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <chrono>



#define LOCAL_SIZE 8

using namespace std;
using namespace cv;
using namespace std::chrono;

typedef cl_float8 float8;
typedef cl_float4 float4;
typedef cl_float3 float3;
typedef cl_float2 float2;

typedef cl_int3 int3;
typedef cl_int2 int2;

typedef cl_uchar8 u8;
typedef cl_uchar3 u3;

struct camera {
	float fx, fy, cx, cy;
};

struct Stereo_Settings_Type1{
	 float3 k1[3], k2[3];
	 float3 R[3], T;
	 float dmin, dmax;
};

struct Stereo_Settings_Type2 {
	float4 p1[3], p2[3];
};

inline int2 make_int2(int x, int y)
{
	int2 tmp;
	tmp.x = x; tmp.y = y;
	return tmp;
}

inline float3 make_float3(float x, float y, float z) 
{
	float3 temp;
	temp.x = x; temp.y = y; temp.z = z;
	return temp;
}


template <typename T, typename u>
inline  T make_vec3(u x, u y, u z)
{
	T temp;
	temp.x = x; temp.y = y; temp.z = z;
	return temp; 
}

#endif
