#include "stdafx.h"
#include "cameraSystem.h"


camera_system::camera_system(float f1, float cx1, float cy1, float f2, float cx2, float cy2, float doff, float dscale, float baseline, int dmin, int dmax)
{
	
	k1 = new float[9];
	k2 = new float[9];

	f1  = f1  * dscale;
	cx1 = cx1 * dscale;
	cy1 = cy1 * dscale;
	
	f2  = f2  * dscale;
	cx2 = cx2 * dscale;
	cy2 = cy2 * dscale;

	doff = doff * dscale;

	// ----------------------
	this->dmin = floor(dmin * dscale);
	this->dmax = ceil(dmax * dscale);

	if (f1 == f2)
	{
		this->zmin = f1 * baseline / (this->dmax + doff);
		this->zmax = ceil((f1 * baseline) / (this->dmin + doff) );
	}

	// --------------Initialize K1
	/**/
	k1[0] = f1;
	k1[1] = 0;
	k1[2] = cx1;

	k1[3] = 0;
	k1[4] = f1;
	k1[5] = cy1;

	k1[6] = 0;
	k1[7] = 0;
	k1[8] = 1;

	// Initialize K2
	k2[0] = f2;
	k2[1] = 0;
	k2[2] = cx2;

	k2[3] = 0;
	k2[4] = f2;
	k2[5] = cy2;

	k2[6] = 0;
	k2[7] = 0;
	k2[8] = 1;

	/**/
	p1 = new float[12];
	p2 = new float[12];

	for (int i = 0 ; i < 3 ; i++)
		for (int j = 0 ; j < 4 ; j++)
		{
			if (i < 3 && j < 3)
			{
				p1[i * 3 + j] = k1[i * 3 + j];
				p2[i * 3 + j] = k2[i * 3 + j];
			}

			if (j == 3)
			{
				p1[i * 3 + j] = 0;
				p2[i * 3 + j] = 0;
			}
		}
	p2[3] = baseline;

	this->bl = baseline;
	this->f    = f1;
	this->cx1  = cx1;
	this->cx2  = cx2;
	this->cy1  = cy1;
	this->cy2  = cy2;
	this->doff = doff;
}


camera_system::~camera_system()
{

}


float* camera_system::getP1()
{
	return p1;
}


float* camera_system::getP2()
{
	return p2;
}
