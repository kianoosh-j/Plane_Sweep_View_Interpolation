#include "stdafx.h"
#include "planeSweep.h"



planeSweep::planeSweep(u3 *img_array, camera_system *settings, int height, int width, int no_views)
{
	this->img_array = img_array;
	this->p1 = settings->getP1();
	this->p2 = settings->getP2();
	this->settings = settings;
	this->height   = height;
	this->width    = width;
	this->no_view  = no_views;
	this->wSize = make_int2(5, 5); // window size is 11 by default

	compile_device_code();
}




planeSweep ::~planeSweep()
{
	/*delete imgZ_dev_out;
	delete imgZ_host_out;
	delete virtual_view_host;
	delete virtual_view_dev;
	*/
	//delete img_array;
}



void planeSweep::compile_device_code()
{
	std::vector<cl::Platform> platforms;
	cl::Platform::get(&platforms);

	auto default_platform = platforms.begin();
	std::vector<cl::Device> devices;
	default_platform->getDevices(CL_DEVICE_TYPE_GPU, &devices);

	cl::Context context(devices);

	std::ifstream kernel_file("clcode.cl");
	std::string src((std::istreambuf_iterator<char>(kernel_file)),
		std::istreambuf_iterator<char>());
	const char* src_new = src.c_str();
	cl::Program program_temp(context, cl::Program::Sources(
		1, std::make_pair(src_new, src.length() + 1)
	));

	this->program = program_temp;

	cl_int err = this->program.build(devices);

	if (err != CL_SUCCESS)
	{
		std::cout << "BUILD ERROR CL. ERROR NO: " << err << std::endl;
		return;
	}
	
	if (err == CL_SUCCESS)
	{
		std::cout << "Build CL SUCCESSFUL" << std::endl;
		return;
	}
		
}




void planeSweep::perform_depth_estimation(int window_size)
{
	if (window_size % 2 == 0)
	{
		std::cout << "ERROR: WINDOW SIZE IS EVEN" << std::endl;
		return;
	}

	if (no_view != 2)
	{
		std::cout << "ERROR: THE SYSTEM IS NOT STEREO" << std::endl;
		return; 
	}
		
	wSize = make_int2((window_size - 1) / 2, (window_size - 1) / 2);

	//auto t1 = high_resolution_clock::now();
	//planeSweepOnHost();
	//auto t2 = high_resolution_clock::now();
	planeSweepOnDevice();
	//auto t3 = high_resolution_clock::now();

	//auto duration_host = duration_cast<milliseconds>(t2 - t1);
	//auto duration_dev = duration_cast<milliseconds>(t3 - t2);


	//cout << "Host time: " << duration_host.count() / 1000.0 << " msec" << std::endl;
	//cout << "Dev time: " << duration_dev.count() / 1000.0 << " msec" << std::endl;

	
	//plane_sweep_compare_host_dev(0, 10, 0, 5);
	//blockMatching();
}


void planeSweep::generate_view(camera *cam_virt, int wSize)
{
	if (wSize % 2 == 0)
	{
		std::cout << "ERROR: WINDOW SIZE IS EVEN" << std::endl;
		return;
	}

	this->wSize = make_int2((wSize - 1) / 2, (wSize - 1) / 2);

	auto t1 = high_resolution_clock::now();
	stereoRenderingOnHost(cam_virt);
	auto t2 = high_resolution_clock::now();
	stereoRenderingOnDevice(cam_virt);
	auto t3 = high_resolution_clock::now();

	auto duration_host = duration_cast<milliseconds>(t2 - t1);
	auto duration_dev = duration_cast<milliseconds>(t3 - t2);
	

	//cout << "Host time: " << duration_host.count() / 1000.0 << " sec" << std::endl;
	//cout << "Total Dev Time: " << duration_dev.count()  << " msec" << std::endl;
	

	stereo_rendering_compare_host_device(0, width, 0, height);
}



// ===============================================================================================
// ======================================== Device Code ==========================================
// ===============================================================================================

void planeSweep::planeSweepOnDevice()
{
	cl::Context context = program.getInfo<CL_PROGRAM_CONTEXT>();
	auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
	cl_int err;

	u3 *L = img_array;
	u3 *R = img_array + height * width;

	cl::Buffer imgL_dev(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		height * width * sizeof(u3), L, &err);

	cl::Buffer imgR_dev(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		height * width * sizeof(u3), R, &err);

	cl::Buffer imgZ_dev(context, CL_MEM_READ_WRITE, height * width * sizeof(float), nullptr, &err);

	if (err != CL_SUCCESS){
		std::cout << "ALLOCATE BUFFER ERROR. ERROR NO: " << err << std::endl;
		return;
	}

	// Kernel
	cl::Kernel kernel(program, "planeSweep_Stereo");
	err = kernel.setArg(0, imgL_dev);
	err = kernel.setArg(1, imgR_dev);
	err = kernel.setArg(2, imgZ_dev);
	err = kernel.setArg(3, this->width);
	err = kernel.setArg(4, this->height);
	err = kernel.setArg(5, settings->getZmin());
	err = kernel.setArg(6, settings->getZmax());
	err = kernel.setArg(7, settings->getFocal());
	err = kernel.setArg(8, settings->getCx1());
	err = kernel.setArg(9, settings->getCy1());
	err = kernel.setArg(10, settings->getCx2());
	err = kernel.setArg(11, settings->getCy2());
	err = kernel.setArg(12, settings->getBaseLine());
	err = kernel.setArg(13, wSize);

	if (err != CL_SUCCESS) {
		std::cout << "SET ARGUMAN ERROR. ERROR NO: " << err << std::endl;
		return;
	}

	// Execute Kernel
	int2 grid_size;
	grid_size.x = ceil((float)width / LOCAL_SIZE) * LOCAL_SIZE;
	grid_size.y = ceil((float)height / LOCAL_SIZE) * LOCAL_SIZE;
	
	cl::CommandQueue queue(context, devices[0]);
	auto t1 = high_resolution_clock::now();
	err = queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(grid_size.x, grid_size.y), cl::NDRange(LOCAL_SIZE, LOCAL_SIZE));
	queue.finish();
	auto t2 = high_resolution_clock::now();
	auto duration_dev = duration_cast<milliseconds>(t2 - t1);
	std::cout << "Dev time: " << duration_dev.count() << " msec" << std::endl;


	// ------------------------------------------------------------------------------------------------
	// Read From Host To Device
	/**/
	imgZ_dev_out = new float[height * width];
	err = queue.enqueueReadBuffer(imgZ_dev, CL_TRUE, 0, height * width * sizeof(float), imgZ_dev_out);
	queue.finish();
	
	if (err != CL_SUCCESS) {
		std::cout << "ERROR: READ FROM BUFFER. ERROR NO: " << err << std::endl;
		return;
	}
	/**/
	//-------------------------------------save the output---------------------------------------------
	// Save the image output:
	/**/
	cv::Mat L_disp(height, width, CV_32FC1);
	cv::Mat L_disp_uc(height, width, CV_8UC1);

	int zMin = settings->getZmin();
	int zMax = settings->getZmax();

	for (int i = 0 ; i < height ; i++)
		for (int j = 0 ; j < width ; j++)
		{
			float z = imgZ_dev_out[i * width + j];
			float z_norm = (z - zMin) / (zMax - zMin);
			uchar z_uchar = (uchar)(z_norm * 255);
			//L_disp.at<float>(i, j) = disp_f;
			L_disp_uc.at<uchar>(i, j) = z_uchar;
		}

	//imwrite("Z_DevV1_W11-delat0point015.jpg", L_disp_uc);
	imwrite("Z_DevV1_W13-delat0point015.jpg", L_disp_uc);
	/**/	
}



void planeSweep::stereoRenderingOnDevice(camera *cam_virt)
{
	cl::Context context = program.getInfo<CL_PROGRAM_CONTEXT>();
	auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
	cl_int err;

	// ------------------------------------------------------------
	u3 *imgL = img_array;
	u3 *imgR = img_array + height * width;
	int imgSZ = height * width;
	
	cl::Buffer imgL_dev(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, imgSZ * sizeof(u3), imgL, &err); 
	cl::Buffer imgR_dev(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, imgSZ * sizeof(u3), imgR, &err);
	cl::Buffer imgVirt_dev(context, CL_MEM_READ_WRITE, imgSZ * sizeof(u3), nullptr, &err);

	if (err != CL_SUCCESS) {
		std::cout << "ALLOCATE BUFFER ERROR. ERROR NO: " << err << std::endl;
		return;
	}
	// ---------------------------------------------------------
	// DEFINE THE KERNEL
	cl::Kernel kernel(program, "imgRendering_Stereo_V2");
	err = kernel.setArg(0, imgL_dev);
	err = kernel.setArg(1, imgR_dev);
	err = kernel.setArg(2, imgVirt_dev);
	err = kernel.setArg(3, height);
	err = kernel.setArg(4, width);
	err = kernel.setArg(5, settings->getFocal());
	err = kernel.setArg(6, settings->getCx1());
	err = kernel.setArg(7, settings->getCy1());
	err = kernel.setArg(8, settings->getCx2());
	err = kernel.setArg(9, settings->getCy2());
	err = kernel.setArg(10, cam_virt->fx);
	err = kernel.setArg(11, cam_virt->cx);
	err = kernel.setArg(12, cam_virt->cy);
	err = kernel.setArg(13, settings->getBaseLine());
	err = kernel.setArg(14, wSize);
	err = kernel.setArg(15, settings->getZmin());
	err = kernel.setArg(16, settings->getZmax());
	err = kernel.setArg(17, 3 * LOCAL_SIZE * 3 * LOCAL_SIZE * sizeof(u3), nullptr);
	err = kernel.setArg(18, 3 * LOCAL_SIZE * 3 * LOCAL_SIZE * sizeof(u3), nullptr);

	
	if (err != CL_SUCCESS)
	{
		std::cout << "ERROR: SET ARGUMAN. ERROR NO: " << err << std::endl;
		return;
	}
	
	//------------------------------------------------------------
	// EXECUTE THE KERNEL
	int2 grid_size;
	grid_size.x = ceil((float)width / LOCAL_SIZE) * LOCAL_SIZE;
	grid_size.y = ceil((float)height / LOCAL_SIZE) * LOCAL_SIZE;

	cl::CommandQueue queue(context, devices[0]);

	auto t1 = high_resolution_clock::now();
	err = queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(grid_size.x, grid_size.y), cl::NDRange(LOCAL_SIZE, LOCAL_SIZE));
	queue.finish();
	auto t2 = high_resolution_clock::now();
	auto duration_dev = duration_cast<milliseconds>(t2 - t1);
	std::cout << "Pure Dev time: " << duration_dev.count() << " msec" << std::endl;

	//-----------------------------------------------------------
	// TESTING
	/**/
	// ----------------------------------------------------------
	// UPDATE FROM HOST TO DEVICE
	virtual_view_dev = new u3[imgSZ];
	err = queue.enqueueReadBuffer(imgVirt_dev, CL_TRUE, 0, imgSZ * sizeof(u3), virtual_view_dev);
	queue.finish();

	if (err != CL_SUCCESS)
	{
		std::cout << "ERROR: READ BUFFER. ERROR NO: " << err << std::endl;
		return;
	}
	// -------------------------------------------------------------
	// SAVE THE OUTPUT IMAGE
	/**/
	cv::Mat img_disk_uc(height, width, CV_8UC3);

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
		{
			u3 pixl = virtual_view_dev[i * width + j];
			img_disk_uc.at<Vec3b>(i, j)[2] = pixl.x;	// R
			img_disk_uc.at<Vec3b>(i, j)[1] = pixl.y;	// G
			img_disk_uc.at<Vec3b>(i, j)[0] = pixl.z;	// B
		}

	imwrite("Intrp_DV4.jpg", img_disk_uc);
	/**/
	
}


//===================================================================================================
//======================================== Host Code ================================================
//===================================================================================================



u3 planeSweep::stereoRenderingPerPixel(camera *cam_virt, int x, int y)
{
	float cost_min = std::numeric_limits<float>::infinity();
	int2 p1_best = make_int2(0, 0);
	int2 p2_best = make_int2(0, 0);


	// Virtual Camera Parameters
	auto f_vrt  = cam_virt->fx;
	auto cx_vrt = cam_virt->cx;
	auto cy_vrt = cam_virt->cy;

	// System Parameters
	auto f   = settings->getFocal();
	auto cx1 = settings->getCx1();
	auto cy1 = settings->getCy1();
	auto cx2 = settings->getCx2();
	auto cy2 = settings->getCy2();
	auto bl  = settings->getBaseLine();

	// -------------------------
	float pos_vrt = 0.45;
	float T1 = bl * pos_vrt;
	float T2 = bl * (1 - pos_vrt);
	//bl = bl / 2;
	
	for (float z = settings->getZmin() ; z <= settings->getZmax() ; z += 0.015)
	{
		float3 xyz_vrt, xyz_left, xyz_right;

		// PROJECT FROM VIRTUAL CAMERA TO THE 3D SPACE
		xyz_vrt.x = (x - cx_vrt) * z / f_vrt;
		xyz_vrt.y = (y - cy_vrt) * z / f_vrt;
		xyz_vrt.z = z;

		// TRANSFORM TO LEFT CAMERA COORDINATE SYSTEM
		xyz_left.x  = xyz_vrt.x + T1;	// Translate to Left
		xyz_left.y  = xyz_vrt.y;
		xyz_left.z  = xyz_vrt.z;

		// TRANSFORM TO RIGHT CAMERA COORDINATE SYSTEM
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
			float cost = SSD_Error(u1, v1, u2, v2);
			if (cost < cost_min)
			{
				cost_min = cost;
				p1_best.x = u1; p1_best.y = v1;
				p2_best.x = u2; p2_best.y = v2;
			}
		}
	}
	
	//debug = p1_best.x;
	/**/
	auto I1 = img_array;
	auto I2 = img_array + height * width;
	auto color_left  = I1[p1_best.y * width + p1_best.x];
	auto color_right = I2[p2_best.y * width + p2_best.x];

	u3 color_best = make_vec3<u3, uchar>((color_left.x + color_right.x) / 2, (color_left.y + color_right.y) / 2,
		(color_left.z + color_right.z) / 2);
	/**/
	//u3 color_best; 
	//color_best = color_right;
	return color_best;
}


void planeSweep::stereoRenderingOnHost(camera *cam_virt)
{
	cv::Mat img_virt_uc(height, width, CV_8UC3);
	
	int zMin = settings->getZmin();
	int zMax = settings->getZmax();
	virtual_view_host = new u3[width * height];

#pragma omp parallel for
	for (int i = 0; i < height ; i++)
		for (int j = 0; j < width ; j++) {
			u3 pixl = stereoRenderingPerPixel(cam_virt, j, i);

			img_virt_uc.at<Vec3b>(i, j)[2] = pixl.x; // R
			img_virt_uc.at<Vec3b>(i, j)[1] = pixl.y; // G
			img_virt_uc.at<Vec3b>(i, j)[0] = pixl.z; // B

			virtual_view_host[i * width + j] = pixl;
		}

	imwrite("intrp_Host.jpg", img_virt_uc);
}


float planeSweep::planeSweepPerPixel(int x, int y)
{

	float cost_min = std::numeric_limits<float>::infinity();
	int d_best = settings->getDmin();	
	float z_best = settings->getZmin();
	

	auto f   = settings->getFocal();
	auto cx1 = settings->getCx1();
	auto cy1 = settings->getCy1();
	auto cx2 = settings->getCx2();
	auto cy2 = settings->getCy2();
	auto bl  = settings->getBaseLine();
	auto doff = cx2 - cx1;

	int uu = -1;

	//for (int d = settings->getDmin() ; d < settings->getDmax() ; d++)
	for (float z = settings->getZmin() ; z <= settings->getZmax() ; z += 0.015)
	{
		//auto z = bl * f / (d + doff);
		float3 xyz;

		// PROJECT FROM LEFT CAMERA TO THE 3D SPACE
		xyz.x = (x - cx1) * z / f;
		xyz.y = (y - cy1) * z / f;
		xyz.z = z;

		// PERFORM SPATIAL TRANSFORM
		xyz.x = xyz.x - bl;	// Translate

		// PROJECT BACK FROM 3D SPACE TO 2D RIGHT CAMERA
		int u = round(xyz.x * f / z + cx2);
		int v = round(xyz.y * f / z + cy2);

		
		if (u > 0)
		{
			float cost = SSD_Error(x, y, u, v);
			if (cost < cost_min)
			{
				cost_min = cost;
				//d_best = d;
				z_best = z;
			}
		}
	}

	return z_best;
}



void planeSweep::planeSweepOnHost()
{
	cv::Mat L_disp(height, width, CV_32FC1);
	cv::Mat L_disp_uc(height, width, CV_8UC1);
	imgZ_host_out = new float[width * height];

	auto zMin = settings->getZmin();
	auto zMax = settings->getZmax();
	auto dMin = settings->getDmin();
	auto dMax = settings->getDmax();

	//#pragma omp parallel for
	for (int i = 0 ; i < height ; i++)
		for (int j = 0 ; j < width ; j++) {

			//int d = planeSweepPerPixel(j, i);
			float z = planeSweepPerPixel(j, i);
			//float z_norm = (z - zMin) / (zMax - zMin);
			//uchar z_uchar = (uchar)((int)(z_norm * 255));
			//L_disp.at<float>(i, j) = disp_f;
			//L_disp_uc.at<uchar>(i, j) = z_uchar;

			imgZ_host_out[i* width + j] = z;

		}

	//imwrite("Z_HostV2_W11-delat0point015.jpg", L_disp_uc);
}




int planeSweep::blockMatchingPerPixel(int x, int y)
{

	float cost_min = std::numeric_limits<float>::infinity();
	int disp_best = settings->getDmin();


	for (int d = settings->getDmin(); d < settings->getDmax(); d++)
	{
		int u = x - d;
		int v = y;
		if (u > 0)
		{
			float cost = SSD_Error(x, y, u, v);
			if (cost < cost_min)
			{
				cost_min = cost;
				disp_best = d;
			}
		}
	}

	return disp_best;
}




void planeSweep::blockMatching()
{
	cv::Mat L_disp(height, width, CV_32FC1);
	cv::Mat L_disp_uc(height, width, CV_8UC1);
	int dMin = settings->getDmin();
	int dMax = settings->getDmax();

#pragma omp parallel for
	for (int i = 0 ; i < height ; i++)
		for (int j = 0 ; j < width ; j++) {

			// CALL BOTH FUNCTION
			int disp = blockMatchingPerPixel(j, i);
			float disp_f = (float)(disp - dMin) / (float)(dMax - dMin);
			uchar disp_uchar = (uchar)((int)(disp_f * 255));
			//L_disp.at<float>(i, j) = disp_f;
			L_disp_uc.at<uchar>(i, j) = disp_uchar;
		}

	imwrite("disp_0.25.jpg", L_disp_uc);
	/**
	namedWindow("My Window", WINDOW_AUTOSIZE);
	imshow("My Window", L_disp_uc);
	waitKey(0);
	/**/
}



float planeSweep::simple_ssd(int x, int y, int u, int v)
{
	u3 *I1 = img_array;
	u3 *I2 = img_array + height * width;

	float ssd = 0.0;

	ssd += float(abs(I1[y * width + x].s0 - I2[v * width + u].s0));
	ssd += float(abs(I1[y * width + x].s1 - I2[v * width + u].s1));
	ssd += float(abs(I1[y * width + x].s2 - I2[v * width + u].s2));

	return ssd;
}


float planeSweep::SSD_Error(int x, int y, int u, int v)
{
	int xstart = x - wSize.x >= 0 ? x - wSize.x : 0;
	int ystart = y - wSize.y >= 0 ? y - wSize.y : 0;

	int xend = x + wSize.x < width  ? x + wSize.x : width  - 1;
	int yend = y + wSize.y < height ? y + wSize.y : height - 1;

	int ustart = u - wSize.x >= 0 ? u - wSize.x : 0;
	int vstart = v - wSize.y >= 0 ? v - wSize.y : 0;

	int uend = u + wSize.x < width  ? u + wSize.x : width - 1;
	int vend = v + wSize.y < height ? v + wSize.y : height - 1;
	
	float ssd = 0;

	u3 *I1 = img_array;
	u3 *I2 = img_array + height * width;


	for (int i = ystart, ii = vstart; (i <= yend) || (ii <= vend); i++, ii++)
	{
		for (int j = xstart, jj = ustart; (j <= xend) || (jj <= uend); j++, jj++)
		{
			ssd += float(abs(I1[i * width + j].s0 - I2[ii * width + jj].s0) );
			ssd += float(abs(I1[i * width + j].s1 - I2[ii * width + jj].s1) );
			ssd += float(abs(I1[i * width + j].s2 - I2[ii * width + jj].s2) );
		}
	}

	return ssd;
}



void planeSweep::plane_sweep_compare_host_dev(int xstart, int xend, int ystart, int yend)
{
	int count = 0;

	for (int i = ystart ; i < yend ; i++)
	{
		for (int j = xstart ; j < xend ; j++)
		{
			int idx = i * width + j;
			auto val_host = imgZ_dev_out[idx];
			auto val_dev  = imgZ_host_out[idx];

			if (val_host != val_dev)
			{
				count++;
				std::cout << "i = " << i << ", j = " << j << ", val_host = " << val_host << ", val_dev = " << val_dev << std::endl;
			}			
		}
	}
}



void planeSweep::stereo_rendering_compare_host_device(int xstart, int xend, int ystart, int yend)
{
	int count = 0;

	for (int i = ystart ; i < yend ; i++)
	{
		for (int j = xstart ; j < xend ; j++)
		{
			int idx = i * width + j;
			u3 val_dev  = virtual_view_dev[idx];
			u3 val_host = virtual_view_host[idx];
			
			if (val_host.x != val_dev.x)
			{
				count++;
				if (count < 20)
					std::cout << "i = " << i << ", j = " << j << ", val_host = " << float(val_host.x) << ", val_dev = " << float(val_dev.x) << std::endl;

			}
			//std::cout << "i = " << i << ", j = " << j << ", val_host = " << float(val_host.x) << ", val_dev = " << float(val_dev.x) << std::endl;
		}
	}

	if (count == 0)
		std::cout << "REDNERING WAS SUCCESSFUL" << std::endl;
	else
		std::cout << "RENDERING ERROR. MISMATCH COUNT: " << count << ", Mismatch Rate: "<<(float)(count)/ (float)(height *width) << std::endl;
}


