#ifndef DEF_H
#define DEF_H
#include<iostream>
#include<opencv2/opencv.hpp>
using namespace std;
using namespace cv;

typedef Mat_<uchar> uMat;
typedef Mat_<Vec3b> uMat3;
typedef Mat_<int>   iMat;
typedef Mat_<float> fMat;
typedef Mat_<Vec2f> fMat2;
typedef Mat_<Vec3f> fMat3;

typedef float f32;
typedef uchar u8;
typedef short s16;
typedef  int   i32;
typedef double d64;


#define EPSILON 0.001
#define MAXITER 300
#define SOR_PARAMETER 1.9
#define GAUSSIAN_SIGMA 0.8


#endif