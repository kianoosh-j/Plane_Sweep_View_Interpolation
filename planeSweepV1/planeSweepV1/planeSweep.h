#include "headers.h"
#include "cameraSystem.h"

//template <typename Settings>
class planeSweep {
public:
	planeSweep(u3 *img_array, camera_system *settings, int height, int width, int no_views);
	~planeSweep();

	void perform_depth_estimation(int window_size);
	void generate_view(camera *cam_virt, int wSize);

	inline float* getimgZ() { return this->imgZ_dev_out; }

private:

	// Variables
	camera_system *settings;
	int height, width, no_view;
	int2 wSize;

	// Arrays
	u3 *img_array;
	float *imgZ_dev_out, *imgZ_host_out;
	u3 *virtual_view_host, *virtual_view_dev;

	// Camera Projection Matrix
	float *p1, *p2;

	//OpenCL
	cl::Program program;
	

	// Methods
	void compile_device_code();

	void planeSweepOnDevice();
	void planeSweepOnHost();
	float planeSweepPerPixel(int x, int y);

	void blockMatching();
	int blockMatchingPerPixel(int x, int y);

	void stereoRenderingOnDevice(camera *cam_virt);
	void stereoRenderingOnHost(camera *cam_virt);
	u3 stereoRenderingPerPixel(camera *cam_virt, int x, int y); 	
	
	float SSD_Error(int x, int y, int u, int v);
	void plane_sweep_compare_host_dev(int xstart, int xend, int ystart, int yend);
	void stereo_rendering_compare_host_device(int xstart, int xend, int ystart, int yend);
	float simple_ssd(int x, int y, int u, int v);
};