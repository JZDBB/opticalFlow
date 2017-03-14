#include <stdio.h>
#include <opencv2\opencv.hpp>
#include <iostream>
using namespace std;
using namespace cv;

//int sum(uchar* p, int mSize);
void padding(uchar** p, int raw, int col);
void translate1(uchar** p, uchar* p0, int col);
void get_grad(uchar** p1, uchar** p2, uchar** u, uchar** v, int raw, int col, int imax);

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
	uchar* current1 = (uchar*)img[0].data;
	uchar* current2 = (uchar*)img[1].data;
	uchar* p10 = (uchar*)malloc(im_size*im_size*sizeof(uchar));
	uchar* p20 = (uchar*)malloc(im_size*im_size*sizeof(uchar));
	uchar** p1 = (uchar**)malloc(im_size*sizeof(uchar*));
	uchar** p2 = (uchar**)malloc(im_size*sizeof(uchar*));

	for (int i = 0; i < im_size*im_size; i++)
	{
		*(p10++) = *(current1++);
	}
	for (int i = 0; i < im_size*im_size; i++)
	{
		*(p20++) = *(current2++);
	}
	p10 = p10 - im_size*im_size;
	p20 = p20 - im_size*im_size;
	//for (int i = 0; i < im_size; i++)
	//{
	//	*(p1 + i) = p10 + im_size*i;
	//}
	//for (int i = 0; i < im_size; i++)
	//{
	//	*(p2 + i) = p20 + im_size*i;
	//}

	translate1(p1, p10, im_size);
	translate1(p2, p20, im_size);

	/*for (int i = 0; i < im_size; i++)
	{
		for (int j = 0; j < im_size; j++)
		{
			if (i < 5 || i>495 || j < 5 || j>495)
			{
				*(*(p1 + i) + j) = 0;
			}
		}
	}
	for (int i = 0; i < im_size; i++)
	{
		for (int j = 0; j < im_size; j++)
		{
			if (i < 5 || i>495 || j < 5 || j>495)
			{
				*(*(p2 + i) + j) = 0;
			}
		}
	}*/

	padding(p1, im_size, im_size);
	padding(p2, im_size, im_size);

	uchar** u = (uchar**)malloc(im_size*sizeof(uchar*));
	uchar** v = (uchar**)malloc(im_size*sizeof(uchar*));
	uchar* u0 = (uchar*)malloc(im_size*im_size*sizeof(uchar));
	uchar* v0 = (uchar*)malloc(im_size*im_size*sizeof(uchar));
	translate1(u, u0, im_size);
	translate1(v, v0, im_size);
	get_grad(p1, p2, u, v, im_size, im_size, 20);
	padding(u, im_size, im_size);
	padding(v, im_size, im_size);

	//for (int i = 0; i < 20; i++)
	//{
	//	for (int j = 0; j < 20; j++)
	//	{
	//		printf("%d,", *(*(u + i) + j));
	//		printf("\n");
	//		printf("%d,", *(*(v + i) + j));
	//	}
	//}




	//Hflow(p1, p2);

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

void get_grad(uchar** p1, uchar** p2, uchar** u, uchar** v, int raw, int col, int imax)
{
	int delta;
	int alpha = 25;
	//uchar* dst0 = (uchar*)malloc(raw*col*sizeof(uchar));
	//uchar* dsx10 = (uchar*)malloc(raw*col*sizeof(uchar));
	//uchar* dsx20 = (uchar*)malloc(raw*col*sizeof(uchar));
	//uchar** dst = (uchar**)malloc(raw*sizeof(uchar*));
	//uchar** dsx1 = (uchar**)malloc(raw*sizeof(uchar*));
	//uchar** dsx2 = (uchar**)malloc(raw*sizeof(uchar*));

	//translate1(dst, dst0, col);
	//translate1(dsx1, dsx10, col);
	//translate1(dsx2, dsx20, col);

	int dst = 0;
	int dsx1 = 0;
	int dsx2 = 0;

	for (int i = 0; i < raw-1; i++)
	{
		for (int j = 0; j < col-1; j++)
		{
			*(*(u + i) + j) = 0;
			*(*(v + i) + j) = 0;
			dst = (int)(*(*(p2 + i + 1) + j + 1) - *(*(p1 + i + 1) + j + 1) + *(*(p2 + i) + j + 1) - *(*(p1 + i) + j + 1) + *(*(p2 + i + 1) + j) - *(*(p1 + i + 1) + j) + *(*(p2 + i) + j) - *(*(p1 + i) + j)) / 4;
			dsx1 = (int)(*(*(p2 + i + 1) + j + 1) - *(*(p2 + i) + j + 1) + *(*(p2 + i +1) + j) - *(*(p2 + i) + j) + *(*(p1 + i + 1) + j + 1) - *(*(p1 + i) + j + 1) + *(*(p1 + i + 1) + j) - *(*(p1 + i) + j)) / 4;
			dsx2 = (int)(*(*(p2 + i + 1) + j + 1) - *(*(p2 + i + 1) + j) + *(*(p2 + i) + j + 1) - *(*(p2 + i) + j) + *(*(p1 + i + 1) + j + 1) - *(*(p1 + i + 1) + j) + *(*(p1 + i) + j + 1) - *(*(p1 + i) + j)) / 4;
			/*if(dsx1!=0)
			{ 
				system("pause");
			}*/
			for (int k = 0; k < imax; k++)
			{
				delta = (int)(dst + *(*(u + i) + j) * dsx1 + *(*(v + i) + j) * dsx2) / (alpha*alpha + dsx1*dsx1 + dsx2*dsx2);
				*(*(u + i) + j) = (uchar)(*(*(u + i) + j) - dsx1 * delta);
				*(*(v + i) + j) = (uchar)(*(*(v + i) + j) - dsx2 * delta);
			}
		}
	}
}

void padding(uchar** p, int raw, int col)
{
	for (int i = 0; i < raw; i++)
	{
		for (int j = 0; j < col; j++)
		{
			if (i < 5 || i>495 || j < 5 || j>495)
			{
				*(*(p + i) + j) = 0;
			}
		}
	}
}

void translate1(uchar** p, uchar* p0, int col)
{
	for (int i = 0; i < col; i++)
	{
		*(p + i) = p0 + col*i;
	}
}