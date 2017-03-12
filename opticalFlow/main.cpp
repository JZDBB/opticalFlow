#include <stdio.h>
#include <opencv2\opencv.hpp>
#include <iostream>
using namespace std;
using namespace cv;
int Hflow(uchar* p1, uchar* p2);
int sum(uchar* p, int mSize);
int main()
{
	//int a = 65;
	//printf("%c\n", a);
	//system("pause");
	int n = 2;
	Mat img[2];
	for (int i = 1; i <= n; i++)
	{
		char filename[10];
		sprintf_s(filename, "%d.bmp", i);
		img[i - 1] = imread(filename, 0);
	}
	//for (int i = 0; i < 4; i++)
	//	for (int j = 0; j < 3; j++)
	//		cout << int(img[0].at<uchar>(i, j)) << endl;
	//system("pause");

	//for (int i = 0; i < 3; i++)
	//{
	//	imshow("", img[i]);
	//	waitkey(0);
	//}
	//return 1;
	int im_size = 500;
	uchar* p1 = (uchar*)malloc(im_size*im_size*sizeof(uchar));
	uchar* p2 = (uchar*)malloc(im_size*im_size*sizeof(uchar));
	//p = (uchar*)malloc(12 * sizeof(uchar));
	uchar* current1 = (uchar*)img[0].data;
	uchar* current2 = (uchar*)img[1].data;
	for (int i = 0; i < im_size*im_size; i++)
	{
		*(p1++) = *(current1++);
	}
	for (int i = 0; i < im_size*im_size; i++)
	{
		*(p2++) = *(current2++);
	}

	//int a = Hflow(p1, p2);

	return 0;
}


//int Hflow(uchar* p1)
//{
//	uchar* p2=NULL;
//	int im_size = 500;
//	uchar Xn[500][500];
//	uchar Xnp1[500][500];
//	for (int i = 0; i < im_size; i++)
//	{
//		for (int j = 0; j < im_size; j++)
//		{
//			Xn[i][j] = *(p1++);
//		}
//	}
//	for (int i = 0; i < im_size; i++)
//	{
//		for (int j = 0; j < im_size; j++)
//		{
//			Xnp1[i][j] = *(p2++);
//		}
//	}
//	
//	uchar dst[500][500], dsx1[500][500], dsx2[500][500], u[500][500], v[500][500];
//	memset(dst, 0, im_size*im_size*sizeof(uchar));
//	memset(dsx1, 0, im_size*im_size*sizeof(uchar));
//	memset(dsx2, 0, im_size*im_size*sizeof(uchar));
//	memset(u, 0, im_size*im_size*sizeof(uchar));
//	memset(v, 0, im_size*im_size*sizeof(uchar));
//	int alpha = 25;
//	int imax = 20;
//	uchar delta;
//	for (int i = 4; i < im_size-5; i++)
//	{
//		for (int j = 4; j < im_size-5; j++)
//		{
//			dst[i][j] = (Xnp1[i + 1][j + 1] - Xn[i + 1][j + 1] + Xnp1[i + 1][j] - Xn[i + 1][j] + Xnp1[i][j + 1] - Xn[i][j + 1] + Xnp1[i][j] - Xn[i][j]) / 4;
//			dsx1[i][j] = (Xnp1[i + 1][j + 1] - Xnp1[i][j + 1] + Xnp1[i + 1][j] - Xnp1[i][j] + Xn[i + 1][j + 1] - Xn[i][j + 1] + Xn[i + 1][j] - Xn[i][j]) / 4;
//			dsx2[i][j] = (Xnp1[i + 1][j + 1] - Xnp1[i + 1][j] + Xnp1[i][j + 1] - Xnp1[i][j] + Xn[i + 1][j + 1] - Xn[i + 1][j] + Xn[i][j + 1] - Xn[i][j]) / 4;
//			for (int k = 0; k < imax; k++)
//			{
//				delta = (dst[i][j] + dsx1[i][j] * u[i][j] + dsx2[i][j] * v[i][j]) / (alpha ^ 2 + dsx1[i][j] ^ 2 + dsx2[i][j] ^ 2);
//				u[i][j] = u[i][j] - delta*dsx1[i][j];
//				v[i][j] = v[i][j] - delta*dsx2[i][j];
//			}
//		}
//	}
//	
//	int im_sizes = 64;
//	int N = 64;
//	uchar us[64][64], vs[64][64], u1[8][8], v1[8][8];
//	memset(us, 0, im_sizes*im_sizes*sizeof(uchar));
//	memset(vs, 0, im_sizes*im_sizes*sizeof(uchar));
//	memset(u1, 0, 64*sizeof(uchar));
//	memset(v1, 0, 64*sizeof(uchar));
//	uchar* pointer1;
//	uchar* pointer2;
//	pointer1 = *u1;
//	pointer2 = *v1;
//	for (int i = 0; i < im_sizes-1; i++)
//	{
//		for (int j = 0; j < 8; j++)
//		{
//			for (int k = 0; k < 8; k++)
//			{
//				u1[j][k] = u[i * 8 + j][i * 8 + k];
//				v1[j][k] = v[i * 8 + j][i * 8 + k];
//				for (int m = 0; m < 64; m++)
//				{
//					for (int n = 0; n < 64; n++)
//					{
//						us[m][n] = (uchar)sum(p1, im_sizes) / N;
//						vs[m][n] = (uchar)sum(p2, im_sizes) / N;
//					}
//				}
//			}
//		}
//	}
//	
//	system("pause");
//	return 0;
//}

//int sum(uchar* p, int mSize)
//{
//	int sum = 0;
//	for (int i = 0; i < mSize; i++)
//	{
//		sum = sum + (int)*(p++);
//	}
//	return sum;
//}