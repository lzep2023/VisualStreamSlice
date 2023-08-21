#include"streamSlice.hpp"

int main(int argc, char* argv[])
{
	int r = 17;
	int c = 33;
	cv::Mat1f x(r, c, CV_32FC1);
	cv::Mat1f y(r, c, CV_32FC1);
	cv::Mat1f u(r, c, CV_32FC1);
	cv::Mat1f v(r, c, CV_32FC1);

	for (int i = 0; i < r; i++)
	for (int j = 0; j < c; j++)
	{
		x(i, j) = j;
		y(i, j) = i;
		u(i, j) = sin(j);
		v(i, j) = cos(i);
	}

	f32 density = 1;
	bool arrows = true;
	cv::Mat3b slImg;
	streamSlice2D(x, y, u, v, density, arrows, slImg);

	cv::imshow("slImg", slImg);
	cv::waitKey(0);
	return 0;
}