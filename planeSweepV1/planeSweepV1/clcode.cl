



float SSD_Error(global uchar3 *L_in, global uchar3 *R_in, int x, int y, int u, int v, int height, int width, int2 wSize);
float simple_ssd(global uchar3 *L_in, global uchar3 *R_in, int x, int y, int u, int v, int height, int width, int2 wSize);
float SSD_Shared(global uchar3 *L_in, global uchar3 *R_in, local uchar3 *L_sh, local uchar3 *R_sh, int x, int y, int u, int v, int height, int width, int2 wSize);


kernel void planeSweep_Stereo(global uchar3 * L_in, global uchar3 * R_in, global float *Z_out, int width, int height, float zMin, float zMax, float f, float cx1, float cy1, float cx2,
 float cy2, float baseline, int2 radious)
{
	
	int x = get_global_id(0);
	int y = get_global_id(1);
	
	if (x >= width || y >= height)
		return;
	
	float cost_min = 9999999999;
	float z_best = zMin;
	float doff = cx2 - cx1;
	
	int uu;
	
	for (float z = zMin ; z <= zMax ; z += 0.015)
	{
		float3 xyz;
		
		// PROJECT FROM LEFT CAMERA TO THE 3D SPACE
		xyz.x = (x - cx1) * z / f;
		xyz.y = (y - cy1) * z / f;
		xyz.z = z; 
		
		// PERFORM SPATIAL TRANSFORM
		xyz.x = xyz.x - baseline;	// Translate

		// PROJECT BACK FROM 3D-SPACE TO 2D RIGHT CAMERA
		int u = round(xyz.x * f / z + cx2);
		int v = round(xyz.y * f / z + cy2);
		
		if (u > 0)
		{
			float cost = SSD_Error(L_in, R_in, x, y, u, v, height, width, radious);
			 if (cost < cost_min)
			{
				cost_min = cost;
				z_best = z;
			}
		}
	}
	
	Z_out[y * width + x] = z_best;
}




kernel void imgRendering_Stereo(global uchar3 *L_in, global uchar3 *R_in, global uchar3 *V_out, int height, int width, float f, float cx1, float cy1, float cx2, float cy2,
float f_vrt, float cx_vrt, float cy_vrt, float bl, int2 radious, float zMin, float zMax)
 {
	int x = get_global_id(0);
	int y = get_global_id(1);
	 
	if (x >= width || y >= height)
		return;
	
	float cost_min = 999999999;
	int2 p1_best = (int2)(0, 0);
	int2 p2_best = (int2)(0, 0);
	
	
	float pos_vrt = 0.95;
	float T1 = bl * pos_vrt;
	float T2 = bl * (1 - pos_vrt);
	
	int debug = 0;
	float costd = 0.0;
	int n = 8;
	for (float z = zMin ; z <= zMax ; z += 0.015)
	{
		float3 xyz_vrt, xyz_left, xyz_right;
		
		// PROJECT FROM VIRTUAL CAMERA TO THE 3D SPACE
		xyz_vrt.x = (x - cx_vrt) * z / f_vrt;
		xyz_vrt.y = (y - cy_vrt) * z / f_vrt;
		xyz_vrt.z = z;
		
		// TRANSFORM TO LEFT CAMERA
		xyz_left.x  = xyz_vrt.x + T1;	// Translate to Left
		xyz_left.y  = xyz_vrt.y;
		xyz_left.z  = xyz_vrt.z;
		
		// TRANSFORM TO RIGHT CAMERA
		xyz_right.x = xyz_vrt.x - T2;	// Translate to Right
		xyz_right.y = xyz_vrt.y;
		xyz_right.z = xyz_vrt.z;
		
		// PROJECT BACK FROM 3D SPACE TO REAL CAMERA 2-D PLANE
		int u1 = round(xyz_left.x * f / z + cx1);
		int v1 = round(xyz_left.y * f / z + cy1);
		
		int u2 = round(xyz_right.x * f / z + cx2);
		int v2 = round(xyz_right.y * f / z + cy2);
		
		
		if (u1 > 0 && u1 < width && u2 > 0 && u2 < width)
		{
			float cost = SSD_Error(L_in, R_in, u1, v1, u2, v2, height, width, radious);
			costd = cost;
			debug = v1;
			if (cost < cost_min)
			{
				cost_min = cost;
				p1_best.x = u1; p1_best.y = v1;
				p2_best.x = u2; p2_best.y = v2;
			}
		}
	}
	
	//debug = p1_best.x;
	uchar3 color_left  = L_in[p1_best.y * width + p1_best.x];
	uchar3 color_right = R_in[p2_best.y * width + p2_best.x];
	uchar3 color_avrg = (uchar3)((color_left.x + color_right.x) / 2, (color_left.y + color_right.y) / 2, (color_left.z + color_right.z) / 2);
	//uchar3 color_avrg = color_left + color_right;
	//color_avrg = color_avrg / 2;
	
	V_out[y * width + x] = color_avrg;
 }


kernel void imgRendering_Stereo_V2(global uchar3 *L_in, global uchar3 *R_in, global uchar3 *V_out, int height, int width, float f, float cx1, float cy1, float cx2, float cy2,
float f_vrt, float cx_vrt, float cy_vrt, float bl, int2 radious, float zMin, float zMax, local uchar3 *L_sh, local uchar3 *R_sh)
{
	int x = get_global_id(0);
	int y = get_global_id(1);
	 
	if (x >= width || y >= height)
		return;
	
	float cost_min = 999999999;
	int2 p1_best = (int2)(0, 0);
	int2 p2_best = (int2)(0, 0);
	
	
	float pos_vrt = 0.45;
	float T1 = bl * pos_vrt;
	float T2 = bl * (1 - pos_vrt);
	
	int debug = 0;
	float costd = 0.0;
	int n = 8;
	for (float z = zMin ; z <= zMax ; z += 0.015)
	{
		float3 xyz_vrt, xyz_left, xyz_right;
		
		// PROJECT FROM VIRTUAL CAMERA TO THE 3D SPACE
		xyz_vrt.x = (x - cx_vrt) * z / f_vrt;
		xyz_vrt.y = (y - cy_vrt) * z / f_vrt;
		xyz_vrt.z = z;
		
		// TRANSFORM TO LEFT CAMERA
		xyz_left.x  = xyz_vrt.x + T1;	// Translate to Left
		xyz_left.y  = xyz_vrt.y;
		xyz_left.z  = xyz_vrt.z;
		
		// TRANSFORM TO RIGHT CAMERA
		xyz_right.x = xyz_vrt.x - T2;	// Translate to Right
		xyz_right.y = xyz_vrt.y;
		xyz_right.z = xyz_vrt.z;
		
		// PROJECT BACK FROM 3D SPACE TO REAL CAMERA 2-D PLANE
		int u1 = round(xyz_left.x * f / z + cx1);
		int v1 = round(xyz_left.y * f / z + cy1);
		
		int u2 = round(xyz_right.x * f / z + cx2);
		int v2 = round(xyz_right.y * f / z + cy2);
		
		
		if (u1 > 0 && u1 < width && u2 > 0 && u2 < width)
		{
			//float cost = SSD_Error(L_in, R_in, u1, v1, u2, v2, height, width, radious);
			float cost = SSD_Shared(L_in, R_in, L_sh, R_sh, u1, v1, u2, v2, height, width, radious);
			costd = cost;
			debug = v1;
			if (cost < cost_min)
			{
				cost_min = cost;
				p1_best.x = u1; p1_best.y = v1;
				p2_best.x = u2; p2_best.y = v2;
			}
		}
	}
	
	//debug = p1_best.x;
	uchar3 color_left  = L_in[p1_best.y * width + p1_best.x];
	uchar3 color_right = R_in[p2_best.y * width + p2_best.x];
	uchar3 color_avrg = (uchar3)((color_left.x + color_right.x) / 2, (color_left.y + color_right.y) / 2, (color_left.z + color_right.z) / 2);
	//uchar3 color_avrg = color_left + color_right;
	//color_avrg = color_avrg / 2;
	//color_avrg = color_right;
	V_out[y * width + x] = color_avrg;
	
}




// ============================================================================================


float ssdTerm(local uchar3 *L, local uchar3 *R, int x, int y)
{
	float ssd = 0.0;
	return ssd;
}


float SSD_Shared(global uchar3 *L_in, global uchar3 *R_in, local uchar3 *L_sh, local uchar3 *R_sh, int x, int y, int u, int v, int height, int width, int2 wSize)
{
	
	int tx = get_local_id(0); 
	int ty = get_local_id(1);
	int cache_width = 3 * get_local_size(0);
	
	int x_start = x - get_local_size(0);
	int y_start = y - get_local_size(1);
	
	int u_start = u - get_local_size(0);
	int v_start = v - get_local_size(1);
	
	
	for (int i = 0 ; i < 9 ; i++)
	{
		int x_offset = (i % 3) * get_local_size(0);
		int y_offset = (i / 3) * get_local_size(0);
		
		int x_idx = x_start + x_offset;
		int y_idx = y_start + y_offset;
		int u_idx = u_start + x_offset;
		int v_idx = v_start + y_offset;
		
		int x_local = tx + x_offset;
		int y_local = ty + y_offset;
		int id_local= y_local * cache_width + x_local;
		
		L_sh[id_local] = (uchar3)(0);
		R_sh[id_local] = (uchar3)(0);
	
		if (x_idx >= 0 && y_idx >= 0 && u_idx >= 0 && v_idx >= 0 && x_idx < width && y_idx < height && u_idx < width && v_idx < height)
		{
			L_sh[id_local] = L_in[y_idx * width + x_idx];
			R_sh[id_local] = R_in[v_idx * width + u_idx];
		}
	
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	
	int tx_start = tx + get_local_size(0);
	int ty_start = ty + get_local_size(1);
	
	
	
	int up = ty_start - wSize.y;
	int indx = up * cache_width + tx_start;
	while ((L_sh[indx].x == 0 && L_sh[indx].y == 0 && L_sh[indx].z == 0) || (R_sh[indx].x == 0 && R_sh[indx].y == 0 && R_sh[indx].z == 0))	
	{
		up++;
		indx = up * cache_width + tx_start;
	}
		
	
	int down = ty_start + wSize.y;
	indx = down * cache_width + tx_start;
	while ((L_sh[indx].x == 0 && L_sh[indx].y == 0 && L_sh[indx].z == 0) || (R_sh[indx].x == 0 && R_sh[indx].y == 0 && R_sh[indx].z == 0))		
	{
		down--;
		indx = down * cache_width + tx_start;
	}
		

	int left = tx_start - wSize.x;
	indx = ty_start * cache_width + left;
	while ((L_sh[indx].x == 0 && L_sh[indx].y == 0 && L_sh[indx].z == 0) || (R_sh[indx].x == 0 && R_sh[indx].y == 0 && R_sh[indx].z == 0))		
	{
		left++;
		indx = ty_start * cache_width + left;
	}
		
	
	int right = tx_start + wSize.x;
	indx = ty_start * cache_width + right;
	while ((L_sh[indx].x == 0 && L_sh[indx].y == 0 && L_sh[indx].z == 0) || (R_sh[indx].x == 0 && R_sh[indx].y == 0 && R_sh[indx].z == 0))		
	{
		right--;
		indx = ty_start * cache_width + right;
	}
		
	
	float ssd = 0;
	for (int i = up ; i <= down ; i++)
		for (int j = left ; j <= right ; j++)
		{
			ssd += (float)(abs(L_sh[i * cache_width + j].s0 - R_sh[i * cache_width + j].s0) );
			ssd += (float)(abs(L_sh[i * cache_width + j].s1 - R_sh[i * cache_width + j].s1) );
			ssd += (float)(abs(L_sh[i * cache_width + j].s2 - R_sh[i * cache_width + j].s2) );
		}
	
	return ssd;
}



float SSD_Error(global uchar3 *L_in, global uchar3 *R_in, int x, int y, int u, int v, int height, int width, int2 wSize)
{
	 float ssd = 0.0;
	
	int xstart = x - wSize.x >= 0 ? x - wSize.x : 0;
	int ystart = y - wSize.y >= 0 ? y - wSize.y : 0;

	int xend = x + wSize.x < width  ? x + wSize.x : width  - 1;
	int yend = y + wSize.y < height ? y + wSize.y : height - 1;

	int ustart = u - wSize.x >= 0 ? u - wSize.x : 0;
	int vstart = v - wSize.y >= 0 ? v - wSize.y : 0;

	int uend = u + wSize.x < width  ? u + wSize.x : width - 1;
	int vend = v + wSize.y < height ? v + wSize.y : height - 1;
	 
	
	for (int i = ystart, ii = vstart; (i <= yend) || (ii <= vend); i++, ii++)
	{
		for (int j = xstart, jj = ustart; (j <= xend) || (jj <= uend); j++, jj++)
		{
			ssd += (float)(abs(L_in[i * width + j].s0 - R_in[ii * width + jj].s0) );
			ssd += (float)(abs(L_in[i * width + j].s1 - R_in[ii * width + jj].s1) );
			ssd += (float)(abs(L_in[i * width + j].s2 - R_in[ii * width + jj].s2) );
		}
	}
	
	return ssd;
}

float simple_ssd(global uchar3 *L_in, global uchar3 *R_in, int x, int y, int u, int v, int height, int width, int2 wSize)
{
	float ssd = 0.0;
	ssd += (float)(abs(L_in[y * width + x].s0 - R_in[v * width + u].s0));
	ssd += (float)(abs(L_in[y * width + x].s1 - R_in[v * width + u].s1));
	ssd += (float)(abs(L_in[y * width + x].s2 - R_in[v * width + u].s2));
	
	return ssd;
}
