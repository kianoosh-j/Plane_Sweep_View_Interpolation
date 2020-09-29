#include "headers.h"

#ifndef CAMERA_SYSTEM
#define CAMERA_SYSTEM

class camera_system {

public: 
	camera_system(float f1, float cx1, float cy1, float f2, float cx2, float cy2, float doff, float dscale, float baseline, int dmin, int dmax);
	~camera_system();

	float* getP1();
	float* getP2();
	inline int getDmin() { return dmin; }
	inline int getDmax() { return dmax; }
	inline float getZmin() { return zmin; }
	inline float getZmax() { return zmax; }
	inline float getFocal() { return f; }
	inline float getBaseLine() { return bl; }
	inline float getCx1() { return cx1; }
	inline float getCx2() { return cx2; }
	inline float getCy1() { return cy1; }
	inline float getCy2() { return cy2; }
	inline float getDoff() { return doff; }
	
	

private: 

	float *k1, *k2;
	float *p1,*p2;
	int dmin, dmax;
	float zmin, zmax;
	float f;
	float bl;
	float cx1, cy1, cx2, cy2;
	float doff;
};
#endif
