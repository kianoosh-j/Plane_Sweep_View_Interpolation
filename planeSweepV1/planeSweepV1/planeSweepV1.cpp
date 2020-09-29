// planeSweepV1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "headers.h"
#include "planeSweep.h"
#include "cameraSystem.h"

void loadImageIn(Mat &in_cam, u3 *in_img, int height, int width)
{
	for (int i = 0 ; i < height ; i++)
		for (int j = 0 ; j < width ; j++)
		{
			in_img[i*width + j].s0 = in_cam.at<Vec3b>(i, j)[2];	// red
			in_img[i*width + j].s1 = in_cam.at<Vec3b>(i, j)[1];	// green
			in_img[i*width + j].s2 = in_cam.at<Vec3b>(i, j)[0]; // blue 
		}
}


void loadImageOut(Mat &out_cam, cl_uchar3 *out_img, int height, int width)
{
	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
		{
			out_cam.at<Vec3b>(i, j)[2] = out_img[i*width + j].s0; // red
			out_cam.at<Vec3b>(i, j)[1] = out_img[i*width + j].s1; // green
			out_cam.at<Vec3b>(i, j)[0] = out_img[i*width + j].s2; // blue
		}
}


bool read_image_files(std::vector<cv::Mat> &img_file, std::string file_name, float scale)
{
	ifstream info_file(file_name);
	string line;

	//namedWindow("My Window", WINDOW_AUTOSIZE); // for test
	
	if (info_file.is_open())
	{
		while (getline(info_file, line))
		{
			auto img_current = imread(line); // Read current image with OpenCV
			cv::Mat img_new;
			auto height_test = img_current.size().height;
			auto width_test = img_current.size().width;
			//cv::Size s(height * scale, width * scale);
			resize(img_current, img_new, cv::Size(), scale, scale);
			img_file.push_back(img_new);
			//imshow("My Window", img_new);
			//waitKey(0);
		}
	}
	else
		return true;

	return false;
}


bool extractRawImageData(u3 *img_array, std::vector<cv::Mat> img_file, int height, int width, int no_view)
{
	for (int i = 0 ; i < no_view ; i++)
	{
		auto img_array_start = img_array + i * height * width;
		loadImageIn(img_file[i], img_array_start, height, width);
	}
	return false;
}



int main()
{

	// ------------------------------------------------------------------------------------------
	// SYSTEM SETTINGS

	int no_view = 2;
	int wSize = 13;
	float baseline = 176.252 * 0.001f; // in meter
	float doffs = 209.059;
	float cx1 = 1445.577f;
	float cy1 = 984.686f;
	float cx2 = 1654.636f;
	float cy2 = 984.686f;
	float f = 4161.221f;
	float dscale = 0.25f;
	float dmin = 33, dmax = 257;	// for the living room
	//float dmin = 34, dmax = 220;
	camera_system *cam_system = new camera_system(f, cx1, cy1, f, cx2, cy2, doffs, dscale, baseline, dmin, dmax);


	// -----------------------------------------------------------------------------------------
	// READ INPUT IMAGES


	std::vector<cv::Mat> img_file;
	if (read_image_files(img_file, "data.txt", dscale))
	{
		std::cerr << "Read file error" << std::endl;
		return 0;
	}

	int height = img_file[0].size().height;
	int width = img_file[0].size().width;

	u3 *img_array = new u3[no_view * height * width];
	
	if (extractRawImageData(img_array, img_file, height, width, no_view))
	{
		std::cerr << "Extract raw image data error" << std::endl;
		return 0;
	}
	


	// ----------------------------------------------------------------------------------------
	

	planeSweep *prob = new planeSweep(img_array, cam_system, height, width, no_view);

	// ESTIMATE THE DEPTH
	//prob->perform_depth_estimation(wSize);

	// GENERATE NEW VIEW 
	camera *cam_virt = new camera;
	cam_virt->fx = cam_system->getFocal();
	cam_virt->fy = cam_system->getFocal();
	cam_virt->cx = cam_system->getCx1();
	cam_virt->cy = cam_system->getCy1();

	prob->generate_view(cam_virt, wSize);
	
	system("pause");
	return 0;
}

