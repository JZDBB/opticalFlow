#include <stdio.h>
#include <opencv2\opencv.hpp>
#include <iostream>
#include <iomanip>  
#include <fstream> 
using namespace std;
using namespace cv;

//int sum(uchar* p, int mSize);
void padding(float** p, int row, int col);
void translate1(float** p, float* p0, int col, int row);
void get_grad(float** p1, float** p2, float** u, float** v, int row, int col, int imax);
float sum_matrix(float** p, int hk, int hl, int wk, int wl, int N);
int round1(float r);

int main()
{
	//int a = 65;
	//printf("%c\n", a);
	//system("pause");
	int n = 2;
	float* img[2];
	int row = 300;
	int col = 400;
	img[0] = (float*)malloc(row * col * sizeof(float));
	img[1] = (float*)malloc(row * col * sizeof(float));
		
	for (int i = 1; i <= n; i++)
	{
		char filename[10];
		sprintf_s(filename, "%d.bmp", i);
		Mat img1 = imread(filename);
		float* current = img[i - 1];
		//bgr
		uchar* p = (uchar*)img1.data;
		for (int j = 0; j < row*col; j++)
		{
			*(current++) = *(p + 3 * j)*double(0.114) + *(p + 3 * j + 1)*double(0.578) + *(p + 3 * j + 2)*double(0.299);
		}
		//验证
		//current = img[i - 1];
		//ofstream outfile("F:\\一宁\\一宁百度同步盘\\DSP\\MATLAB\\data1.txt", ofstream::out);
		//for (int k = 0; k < row; k++)
		//{
		//	for (int j = 0; j < col; j++)
		//	{
		//		outfile << (float)*(current++) << ",";
		//	}
		//	outfile << "\n";
		//}
		//outfile.close();
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
	
	//float* current1 = (float*)img[0];
	//float* current2 = (float*)img[1];
	float* p10 = img[0];
	float* p20 = img[1];
	float** p1 = (float**)malloc(row*sizeof(float*));
	float** p2 = (float**)malloc(row*sizeof(float*));

	//for (int i = 0; i < im_size*im_size; i++)
	//{
	//	*(p10++) = (float)*(current1++);
	//}
	//for (int i = 0; i < im_size*im_size; i++)
	//{
	//	*(p20++) = (float)*(current2++);
	//}
	//p10 = p10 - im_size*im_size;
	//p20 = p20 - im_size*im_size;
	//for (int i = 0; i < im_size; i++)
	//{
	//	*(p1 + i) = p10 + im_size*i;
	//}
	//for (int i = 0; i < im_size; i++)
	//{
	//	*(p2 + i) = p20 + im_size*i;
	//}

	translate1(p1, p10, col, row);
	translate1(p2, p20, col, row);

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

	//padding(p1, im_size, im_size);
	//padding(p2, im_size, im_size);

	float** u = (float**)malloc(row * sizeof(float*));
	float** v = (float**)malloc(row * sizeof(float*));
	float* u_data = (float*)malloc(row * col * sizeof(float));
	float* v_data = (float*)malloc(row * col * sizeof(float));
	translate1(u, u_data, col, row);
	translate1(v, v_data, col, row);
	get_grad(p1, p2, u, v, row, col, 20);

	//for (int i = 0; i < 20; i++)
	//{
	//	for (int j = 0; j < 20; j++)
	//	{
	//		printf("%d,", *(*(u + i) + j));
	//		printf("\n");
	//		printf("%d,", *(*(v + i) + j));
	//	}
	//}

	//验证
	ofstream outfile("F:\\一宁\\一宁百度同步盘\\DSP\\MATLAB\\out1.txt",ofstream::out);
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			outfile << (float)*(*(u + i) + j) << ",";
		}
		outfile << "\n";
	}
	outfile.close();
	
	outfile.open("F:\\一宁\\一宁百度同步盘\\DSP\\MATLAB\\out2.txt", ofstream::out);
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			outfile << (float)*(*(v + i) + j) << ",";
		}
		outfile << "\n";
	}
	outfile.close();

	int xskip = round1((float)row/64);
	int row1 = row / xskip;
	int col1 = col / xskip;
	int N = xskip*xskip;

	float** us = (float**)malloc(row1 * sizeof(float*));
	float** vs = (float**)malloc(row1 * sizeof(float*));
	float* us_data = (float*)malloc(row1 * col1 * sizeof(float));
	float* vs_data = (float*)malloc(row1 * col1 * sizeof(float));
	translate1(us, us_data, col1, row1);
	translate1(vs, vs_data, col1, row1);

	for (int i = 0; i < row1-1; i++)
	{
		for (int j = 0; j < col1-1; j++)
		{
			int hk = i*xskip + 1;
			int hl = i*xskip + xskip;
			int wk = j*xskip + 1;
			int wl = j*xskip + xskip;
			*(*(us + i) + j) = sum_matrix(u, hk, hl, wk, wl, N);
			*(*(vs + i) + j) = sum_matrix(v, hk, hl, wk, wl, N);
		}
	}

	//验证
	outfile.open("F:\\一宁\\一宁百度同步盘\\DSP\\MATLAB\\outs1.txt",ofstream::out);
	for (int i = 0; i < row1; i++)
	{
		for (int j = 0; j < col1; j++)
		{
			outfile << (float)*(*(us + i) + j) << ",";
		}
		outfile << "\n";
	}
	outfile.close();
	
	outfile.open("F:\\一宁\\一宁百度同步盘\\DSP\\MATLAB\\outs2.txt", ofstream::out);
	for (int i = 0; i < row1; i++)
	{
		for (int j = 0; j < col1; j++)
		{
			outfile << (float)*(*(vs + i) + j) << ",";
		}
		outfile << "\n";
	}
	outfile.close();

	return 0;
}


//int Hflow(uchar* p1)//栈溢出弃用，数组存在栈
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

void get_grad(float** p1, float** p2, float** u, float** v, int row, int col, int imax)
{
	float delta;
	int alpha = 625;
	//uchar* dst0 = (uchar*)malloc(raw*col*sizeof(uchar));
	//uchar* dsx10 = (uchar*)malloc(raw*col*sizeof(uchar));
	//uchar* dsx20 = (uchar*)malloc(raw*col*sizeof(uchar));
	//uchar** dst = (uchar**)malloc(raw*sizeof(uchar*));
	//uchar** dsx1 = (uchar**)malloc(raw*sizeof(uchar*));
	//uchar** dsx2 = (uchar**)malloc(raw*sizeof(uchar*));

	//translate1(dst, dst0, col);
	//translate1(dsx1, dsx10, col);
	//translate1(dsx2, dsx20, col);

	float dst = 0;
	float dsx1 = 0;
	float dsx2 = 0;

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			*(*(u + i) + j) = 0;
			*(*(v + i) + j) = 0;
			if (i < 4 || i>row-5 || j < 4 || j>col-5)
			{
				
				dst = 0;
				dsx1 = 0;
				dsx2 = 0;
			}
			else {
				
				dst = (float)(*(*(p2 + i + 1) + j + 1) - *(*(p1 + i + 1) + j + 1) + *(*(p2 + i) + j + 1) - *(*(p1 + i) + j + 1) + *(*(p2 + i + 1) + j) - *(*(p1 + i + 1) + j) + *(*(p2 + i) + j) - *(*(p1 + i) + j)) / 4;
				dsx2 = (float)(*(*(p2 + i + 1) + j + 1) - *(*(p2 + i) + j + 1) + *(*(p2 + i + 1) + j) - *(*(p2 + i) + j) + *(*(p1 + i + 1) + j + 1) - *(*(p1 + i) + j + 1) + *(*(p1 + i + 1) + j) - *(*(p1 + i) + j)) / 4;
				dsx1 = (float)(*(*(p2 + i + 1) + j + 1) - *(*(p2 + i + 1) + j) + *(*(p2 + i) + j + 1) - *(*(p2 + i) + j) + *(*(p1 + i + 1) + j + 1) - *(*(p1 + i + 1) + j) + *(*(p1 + i) + j + 1) - *(*(p1 + i) + j)) / 4;

			}
			
			for (int k = 0; k < imax; k++)
			{
				delta = (dst + *(*(u + i) + j) * dsx1 + *(*(v + i) + j) * dsx2) / (alpha + dsx1*dsx1 + dsx2*dsx2);
				*(*(u + i) + j) = (float)(*(*(u + i) + j) - dsx1 * delta);
				*(*(v + i) + j) = (float)(*(*(v + i) + j) - dsx2 * delta);
				//if (*(*(u + i) + j) != 0) {
				//	getchar();
				//}
			}
		}
	}
}

//
//void padding(float** p, int row, int col)
//{
//	for (int i = 0; i < row; i++)
//	{
//		for (int j = 0; j < col; j++)
//		{
//			if (i < 5 || i>495 || j < 5 || j>495)
//			{
//				*(*(p + i) + j) = 0;
//			}
//		}
//	}
//}

void translate1(float** p, float* p0, int col, int row)
{
	for (int i = 0; i < row; i++)
	{
		*(p + i) = p0 + col*i;
	}
}

int round1(float r)
{
	return (int)(r - int(r)) > 0.5 ? floor(r) : ceil(r);
}

float sum_matrix(float** p, int hk, int hl, int wk, int wl, int N)
{
	float sum = 0;
	for (int i = hk; i < hl; i++)
	{
		for (int j = wk; j < wl; j++)
		{
			sum = sum + *(*(p + i) + j) / N;
		}
	}
	return sum;
}