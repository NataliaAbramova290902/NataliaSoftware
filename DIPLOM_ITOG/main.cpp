//===============================================================================
//
// Microsoft Press
// C++ AMP: Accelerated Massive Parallelism with Microsoft Visual C++
//
//===============================================================================
// Copyright (c) 2012-2013 Ade Miller & Kate Gregory.  All rights reserved.
// This code released under the terms of the 
// Microsoft Public License (Ms-PL), http://ampbook.codeplex.com/license.
//
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//===============================================================================

// Demonstrates how to use C++ AMP to do matrix multiply.


#include <iostream>
#include <algorithm>
#include <random>
#include <amp.h>
#include <amp_math.h>
#include <fstream>
#include <clocale>
#include <ppl.h>
#include <vector>
#include <tuple>
#include <algorithm>
#include <concurrent_vector.h>
#include <ctime>
#include "Timer.h"


using namespace concurrency;
using namespace std;

//--------------------------------------------------------------------------------------
//  GPU simple matrix multiply
//--------------------------------------------------------------------------------------

struct Point 
{
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

//TODO: ��������� ��� ��������� ����� �����������, ���� �� �� ����� ���� ������ ����
struct TwoPoints
{
	double onePointZ;
	double twoPointZ;	
};


struct Line
{
	double x = 0;
	double y = 0;
};

//��������� ��� �������� ����� ������
struct Capsule
{
	Point onePoint;
	Point twoPoint;
	
};


void GPUCalc_3(int StepX, int StepY, double R, vector<Capsule>& ScopeCenter, vector<TwoPoints>& NeedIntersections)
{
	
	array_view<const Capsule, 1 > SC(ScopeCenter.size(), ScopeCenter);
	array_view<TwoPoints, 2> NI((StepX + 1) , (StepY + 1), NeedIntersections);
	
	for (int s = 0; s < ScopeCenter.size(); s++)
	{

		//������� ������� ���������	

		//������ ������� ������ ���� OX
		Point Vector;
		Vector.x = SC[s].twoPoint.x - SC[s].onePoint.x;
		Vector.y = SC[s].twoPoint.y - SC[s].onePoint.y;
		Vector.z = SC[s].twoPoint.z - SC[s].onePoint.z;

		double SQ_Vector = abs(sqrt((Vector.y * Vector.y) + (Vector.x * Vector.x) + (Vector.z * Vector.z)));

		Point Vector_n;
		Vector_n.x = Vector.x / SQ_Vector;
		Vector_n.y = Vector.y / SQ_Vector;
		Vector_n.z = Vector.z / SQ_Vector;
		//������� ���� ��������
		double SQ_Vector_n = abs(sqrt((Vector_n.y * Vector_n.y) + (Vector_n.x * Vector_n.x)));

		double SinA;
		double CosA;

		if (SQ_Vector_n == 0)
		{
			SinA = 0;
			CosA = 0;
		}
		else
		{
			SinA = Vector_n.y / SQ_Vector_n;
			CosA = Vector_n.x / SQ_Vector_n;
		}
		
		double SinB = Vector_n.z;
		double CosB = SQ_Vector_n;

		//����� ����� �������
		Point SC_One;
		SC_One.x = 0;
		SC_One.y = 0;
		SC_One.z = 0;

		Point SC_Two;
		SC_Two.x = abs(sqrt((Vector.x * Vector.x) + (Vector.y * Vector.y) + (Vector.z * Vector.z)));
		SC_Two.y = 0;
		SC_Two.z = 0;

		//���������� �������������
		int NI1y;
		int NI2y;
		int NI1x;
		int NI2x;

		if (SC[s].onePoint.y < SC[s].twoPoint.y)
		{
			NI2y = (fast_math::ceil(SC[s].twoPoint.y + R) + 1);
			NI1y = (fast_math::floor(SC[s].onePoint.y - R));
		}
		else
		{
			NI2y = (fast_math::ceil(SC[s].onePoint.y + R) + 1);
			NI1y = (fast_math::floor(SC[s].twoPoint.y - R));
		}

		if (SC[s].onePoint.x < SC[s].twoPoint.x)
		{
			NI2x = (fast_math::ceil(SC[s].twoPoint.x + R) + 1);
			NI1x = (fast_math::floor(SC[s].onePoint.x - R));
		}

		else
		{
			NI2x = (fast_math::ceil(SC[s].onePoint.x + R) + 1);
			NI1x = (fast_math::floor(SC[s].twoPoint.x - R));
		}


		//��� ������� ������� ��������� � �������
		concurrency::extent<2> e(NI2y - NI1y, NI2x - NI1x);
		parallel_for_each(e, [=](index<2> idx) restrict(amp)
		{
			//�������� ����������
			int j = (NI1x + idx[1]); 
			int i = (NI1y + idx[0]);

			//Point NI_New1;
			double NI_New1_x = static_cast<double> (j);
			double NI_New1_y = static_cast<double> (i);
			double NI_New1_z = -100000 * R;

			//Point NI_New2;
			double NI_New2_x = static_cast<double> (j);
			double NI_New2_y = static_cast<double> (i);
			double NI_New2_z = 100000 * R;

			//���������� ����� ��������� �������
			//������ �����
			
			//������� � ����� ������ ���������
			NI_New1_x = NI_New1_x - SC[s].onePoint.x;
			NI_New1_y = NI_New1_y - SC[s].onePoint.y;
			NI_New1_z = NI_New1_z - SC[s].onePoint.z;

			//������� ������ ��� OZ
			NI_New1_x = (NI_New1_x * CosA) + (NI_New1_y * (-SinA)) + (NI_New1_z * 0);
			NI_New1_y = (NI_New1_x * SinA) + (NI_New1_y * CosA) + (NI_New1_z * 0);
			NI_New1_z = (NI_New1_x * 0) + (NI_New1_y * 0) + (NI_New1_z * 1);

			//������� ������ ��� OY
			NI_New1_x = (NI_New1_x * CosB) + (NI_New1_y * 0) + (NI_New1_z * (-SinB));
			NI_New1_y = (NI_New1_x * 0) + (NI_New1_y * 1) + (NI_New1_z * 0);
			NI_New1_z = (NI_New1_x * SinB) + (NI_New1_y * 0) + (NI_New1_z * CosB);

			//������ �����
			//������� � ����� ������ ���������
			NI_New2_x = NI_New2_x - SC[s].onePoint.x;
			NI_New2_y = NI_New2_y - SC[s].onePoint.y;
			NI_New2_z = NI_New2_z - SC[s].onePoint.z;

			//������� ������ ��� OZ
			NI_New2_x = (NI_New2_x * CosA) + (NI_New2_y * (-SinA)) + (NI_New2_z * 0);
			NI_New2_y = (NI_New2_x * SinA) + (NI_New2_y * CosA) + (NI_New2_z * 0);
			NI_New2_z = (NI_New2_x * 0) + (NI_New2_y * 0) + (NI_New2_z * 1);

			//������� ������ ��� OY
			NI_New2_x = (NI_New2_x * CosB) + (NI_New2_y * 0) + (NI_New2_z * (-SinB));
			NI_New2_y = (NI_New2_x * 0) + (NI_New2_y * 1) + (NI_New2_z * 0);
			NI_New2_z = (NI_New2_x * SinB) + (NI_New2_y * 0) + (NI_New2_z * CosB);
				
			double A;
			double B;
			double C;
			double D;

			double NI_New1_x_n;
			double NI_New2_x_n;
			double NI_New1_y_n;
			double NI_New2_y_n;
			double NI_New1_z_n;
			double NI_New2_z_n;


			//������ �����
			 {
				A = ((NI_New2_x - NI_New1_x) * (NI_New2_x - NI_New1_x)) + ((NI_New2_y - NI_New1_y)*(NI_New2_y - NI_New1_y)) + ((NI_New2_z - NI_New1_z) * (NI_New2_z - NI_New1_z));
				B = 2 * (((NI_New2_x - NI_New1_x) * (NI_New1_x - SC_One.x)) + ((NI_New2_y - NI_New1_y) * (NI_New1_y - SC_One.y)) + ((NI_New2_z - NI_New1_z) * (NI_New1_z - SC_One.z)));
				C = ((NI_New1_x - SC_One.x) * (NI_New1_x - SC_One.x)) + ((NI_New1_y - SC_One.y) * (NI_New1_y - SC_One.y)) + ((NI_New1_z - SC_One.z) * (NI_New1_z - SC_One.z)) - (R * R);

				D = (B * B) - 4 * A * C;

				if (D > 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					double t2 = (-B - SQ) / (2 * A);

					//������ ����� �����������
					NI_New1_x_n = NI_New1_x + t1 * (NI_New2_x - NI_New1_x);
					NI_New1_y_n = NI_New1_y + t1 * (NI_New2_y - NI_New1_y);
					NI_New1_z_n = NI_New1_z + t1 * (NI_New2_z - NI_New1_z);

					//������ ����� �����������
					NI_New2_x_n = NI_New1_x + t2 * (NI_New2_x - NI_New1_x);
					NI_New2_y_n = NI_New1_y + t2 * (NI_New2_y - NI_New1_y);
					NI_New2_z_n = NI_New1_z + t2 * (NI_New2_z - NI_New1_z);

					//������������� ���������� �������

					//������ �����

					//������� ������ ��� OY
					NI_New1_z_n = (NI_New1_x_n * (-SinB)) + (NI_New1_y_n * 0) + (NI_New1_z_n * CosB);

					//������� ������ ��� OZ
					NI_New1_z_n = (NI_New1_x_n * 0) + (NI_New1_y_n * 0) + (NI_New1_z_n * 1);

					//������� � ����� ������ ���������
					NI_New1_z_n = NI_New1_z_n + SC[s].onePoint.z;

					//������ �����
					//������� ������ ��� OY
					NI_New2_z_n = (NI_New2_x_n * (-SinB)) + (NI_New2_y_n * 0) + (NI_New2_z_n * CosB);

					//������� ������ ��� OZ
					NI_New2_z_n = (NI_New2_x_n * 0) + (NI_New2_y_n * 0) + (NI_New2_z_n * 1);

					//������� � ����� ������ ���������
					NI_New2_z_n = NI_New2_z_n + SC[s].onePoint.z;

					if (NI_New1_z_n < NI_New2_z_n)
					{
						if (NI[i][j].onePointZ > NI_New1_z_n)
							NI[i][j].onePointZ = NI_New1_z_n;

						if (NI[i][j].twoPointZ < NI_New2_z_n)
							NI[i][j].twoPointZ = NI_New2_z_n;
					}
					else
					{
						if (NI[i][j].onePointZ > NI_New2_z_n)
							NI[i][j].onePointZ = NI_New2_z_n;

						if (NI[i][j].twoPointZ < NI_New1_z_n)
							NI[i][j].twoPointZ = NI_New1_z_n;
					}
				}

				else if (D == 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					//������ ����� �����������
					NI_New1_x_n = NI_New1_x + t1 * (NI_New2_x - NI_New1_x);
					NI_New1_y_n = NI_New1_y + t1 * (NI_New2_y - NI_New1_y);
					NI_New1_z_n = NI_New1_z + t1 * (NI_New2_z - NI_New1_z);

					//������ ����� �����������
					NI_New2_x_n = NI_New1_x_n;
					NI_New2_y_n = NI_New1_y_n;
					NI_New2_z_n = NI_New1_z_n;


					//������������� ���������� �������

					//������ �����

					//������� ������ ��� OY
					NI_New1_z_n = (NI_New1_x_n * (-SinB)) + (NI_New1_y_n * 0) + (NI_New1_z_n * CosB);

					//������� ������ ��� OZ
					NI_New1_z_n = (NI_New1_x_n * 0) + (NI_New1_y_n * 0) + (NI_New1_z_n * 1);

					//������� � ����� ������ ���������
					NI_New1_z_n = NI_New1_z_n + SC[s].onePoint.z;

					if (NI[i][j].onePointZ > NI_New1_z_n)
						NI[i][j].onePointZ = NI_New1_z_n;

					if (NI[i][j].twoPointZ < NI_New1_z_n)
						NI[i][j].twoPointZ = NI_New1_z_n;
				}
			}
			//������ �����
			 {
				A = ((NI_New2_x - NI_New1_x) * (NI_New2_x - NI_New1_x)) + ((NI_New2_y - NI_New1_y) * (NI_New2_y - NI_New1_y)) + ((NI_New2_z - NI_New1_z) * (NI_New2_z - NI_New1_z));
				B = 2 * (((NI_New2_x - NI_New1_x) * (NI_New1_x - SC_Two.x)) + ((NI_New2_y - NI_New1_y) * (NI_New1_y - SC_Two.y)) + ((NI_New2_z - NI_New1_z) * (NI_New1_z - SC_Two.z)));
				C = ((NI_New1_x - SC_Two.x) * (NI_New1_x - SC_Two.x)) + ((NI_New1_y - SC_Two.y) * (NI_New1_y - SC_Two.y)) + ((NI_New1_z - SC_Two.z) * (NI_New1_z - SC_Two.z)) - (R * R);

				D = (B * B) - 4 * A * C;

				if (D > 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					double t2 = (-B - SQ) / (2 * A);

					//������ ����� �����������
					NI_New1_x_n = NI_New1_x + t1 * (NI_New2_x - NI_New1_x);
					NI_New1_y_n = NI_New1_y + t1 * (NI_New2_y - NI_New1_y);
					NI_New1_z_n = NI_New1_z + t1 * (NI_New2_z - NI_New1_z);

					//������ ����� �����������
					NI_New2_x_n = NI_New1_x + t2 * (NI_New2_x - NI_New1_x);
					NI_New2_y_n = NI_New1_y + t2 * (NI_New2_y - NI_New1_y);
					NI_New2_z_n = NI_New1_z + t2 * (NI_New2_z - NI_New1_z);

					//������������� ���������� �������

					//������ �����

					//������� ������ ��� OY
					NI_New1_z_n = (NI_New1_x_n * (-SinB)) + (NI_New1_y_n * 0) + (NI_New1_z_n * CosB);

					//������� ������ ��� OZ
					NI_New1_z_n = (NI_New1_x_n * 0) + (NI_New1_y_n * 0) + (NI_New1_z_n * 1);

					//������� � ����� ������ ���������
					NI_New1_z_n = NI_New1_z_n + SC[s].onePoint.z;

					//������ �����
					//������� ������ ��� OY
					NI_New2_z_n = (NI_New2_x_n * (-SinB)) + (NI_New2_y_n * 0) + (NI_New2_z_n * CosB);

					//������� ������ ��� OZ
					NI_New2_z_n = (NI_New2_x_n * 0) + (NI_New2_y_n * 0) + (NI_New2_z_n * 1);

					//������� � ����� ������ ���������
					NI_New2_z_n = NI_New2_z_n + SC[s].onePoint.z;

					if (NI_New1_z_n < NI_New2_z_n)
					{
						if (NI[i][j].onePointZ > NI_New1_z_n)
							NI[i][j].onePointZ = NI_New1_z_n;

						if (NI[i][j].twoPointZ < NI_New2_z_n)
							NI[i][j].twoPointZ = NI_New2_z_n;
					}
					else
					{
						if (NI[i][j].onePointZ > NI_New2_z_n)
							NI[i][j].onePointZ = NI_New2_z_n;

						if (NI[i][j].twoPointZ < NI_New1_z_n)
							NI[i][j].twoPointZ = NI_New1_z_n;
					}
				}

				else if (D == 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);

					//������ ����� �����������
					NI_New1_x_n = NI_New1_x + t1 * (NI_New2_x - NI_New1_x);
					NI_New1_y_n = NI_New1_y + t1 * (NI_New2_y - NI_New1_y);
					NI_New1_z_n = NI_New1_z + t1 * (NI_New2_z - NI_New1_z);


					//������������� ���������� �������

					//������ �����

					//������� ������ ��� OY
					NI_New1_z_n = (NI_New1_x_n * (-SinB)) + (NI_New1_y_n * 0) + (NI_New1_z_n * CosB);

					//������� ������ ��� OZ
					NI_New1_z_n = (NI_New1_x_n * 0) + (NI_New1_y_n * 0) + (NI_New1_z_n * 1);

					//������� � ����� ������ ���������
					NI_New1_z_n = NI_New1_z_n + SC[s].onePoint.z;

					if (NI[i][j].onePointZ > NI_New1_z_n)
						NI[i][j].onePointZ = NI_New1_z_n;

					if (NI[i][j].twoPointZ < NI_New1_z_n)
						NI[i][j].twoPointZ = NI_New1_z_n;
				}
			}

			//�������
			  {
				A = ((NI_New2_y - NI_New1_y) * (NI_New2_y - NI_New1_y)) + ((NI_New2_z - NI_New1_z) * (NI_New2_z - NI_New1_z));
				B = 2 * (((NI_New2_y - NI_New1_y) * (NI_New1_y - SC_One.y)) + ((NI_New2_z - NI_New1_z) * (NI_New1_z - SC_One.z)));
				C = ((NI_New1_y - SC_One.y) * (NI_New1_y - SC_One.y)) + ((NI_New1_z - SC_One.z) * (NI_New1_z - SC_One.z)) - (R * R);

				D = (B * B) - 4 * A * C;

				if (D > 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					double t2 = (-B - SQ) / (2 * A);

					//������ ����� �����������
					NI_New1_x_n = NI_New1_x + t1 * (NI_New2_x - NI_New1_x);
					NI_New1_y_n = NI_New1_y + t1 * (NI_New2_y - NI_New1_y);
					NI_New1_z_n = NI_New1_z + t1 * (NI_New2_z - NI_New1_z);

					//������ ����� �����������
					NI_New2_x_n = NI_New1_x + t2 * (NI_New2_x - NI_New1_x);
					NI_New2_y_n = NI_New1_y + t2 * (NI_New2_y - NI_New1_y);
					NI_New2_z_n = NI_New1_z + t2 * (NI_New2_z - NI_New1_z);


					//���������� ����� �������� �� ������
					int index_c = 0; // ���������� ����� ����� �� �����: 0 - ��� ����� ���������� �������, 1 - ������ �� ������, 2 - ������ �� ������, 3 - ��� �� ������

					if (NI_New1_x_n < SC_One.x || NI_New1_x_n > SC_Two.x)
					{
						index_c = index_c + 1;
					}

					if (NI_New2_x_n < SC_One.x || NI_New2_x_n > SC_Two.x)
					{
						index_c = index_c + 2;
					}

					if (index_c == 0 || index_c == 2) //���� ������ ������
					{
						//������������� ���������� �������

						//������ �����

						//������� ������ ��� OY
						NI_New1_z_n = NI_New1_x_n * (-SinB) + NI_New1_y_n * 0 + NI_New1_z_n * CosB;

						//������� ������ ��� OZ
						NI_New1_z_n = NI_New1_x_n * 0 + NI_New1_y_n * 0 + NI_New1_z_n * 1;

						//������� � ����� ������ ���������
						NI_New1_z_n = NI_New1_z_n + SC[s].onePoint.z;
					}


					if (index_c == 0 || index_c == 1) // ���� ������ ������
					{
						//������ �����
						//������� ������ ��� OY
						NI_New2_z_n = NI_New2_x_n * (-SinB) + NI_New2_y_n * 0 + NI_New2_z_n * CosB;

						//������� ������ ��� OZ
						NI_New2_z_n = NI_New2_x_n * 0 + NI_New2_y_n * 0 + NI_New2_z_n * 1;

						//������� � ����� ������ ���������
						NI_New2_z_n = NI_New2_z_n + SC[s].onePoint.z;
					}

					if (index_c == 1) //��������� ������ ���� ��� �� ������
					{
						NI_New1_x_n = NI_New2_x_n;
						NI_New1_y_n = NI_New2_y_n;
						NI_New1_z_n = NI_New2_z_n;
					}
					else if (index_c == 2) // ��������� ������ ���� ��� �� ������
					{
						NI_New2_x_n = NI_New1_x_n;
						NI_New2_y_n = NI_New1_y_n;
						NI_New2_z_n = NI_New1_z_n;
					}
					else if (index_c == 3)// ��������� ���
					{
						NI_New1_z_n = NI[i][j].onePointZ;
						NI_New2_z_n = NI[i][j].twoPointZ;
					}

					if (NI_New1_z_n != 100000 * R &&  NI_New2_z_n != -100000 * R)
					{
						if (NI_New1_z_n < NI_New2_z_n)
							{
								if (NI[i][j].onePointZ > NI_New1_z_n)
									NI[i][j].onePointZ = NI_New1_z_n;

								if (NI[i][j].twoPointZ < NI_New2_z_n)
									NI[i][j].twoPointZ = NI_New2_z_n;
							}
						else
						{
							if (NI[i][j].onePointZ > NI_New2_z_n)
								NI[i][j].onePointZ = NI_New2_z_n;

							if (NI[i][j].twoPointZ < NI_New1_z_n)
								NI[i][j].twoPointZ = NI_New1_z_n;
						}
					}
				}

				else if (D == 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);

					//������ ����� �����������
					NI_New1_x_n = NI_New1_x + t1 * (NI_New2_x - NI_New1_x);
					NI_New1_y_n = NI_New1_y + t1 * (NI_New2_y - NI_New1_y);
					NI_New1_z_n = NI_New1_z + t1 * (NI_New2_z - NI_New1_z);

					//������ ����� �����������
					NI_New2_x_n = NI_New1_x_n;
					NI_New2_y_n = NI_New1_y_n;
					NI_New2_z_n = NI_New1_z_n;

					if (NI_New1_x_n >= SC_One.x || NI_New1_x_n <= SC_Two.x)
					{
						//������������� ���������� �������

						//������ �����

						//������� ������ ��� OY
						NI_New1_z_n = NI_New1_x_n * (-SinB) + NI_New1_y_n * 0 + NI_New1_z_n * CosB;

						//������� ������ ��� OZ
						NI_New1_z_n = NI_New1_x_n * 0 + NI_New1_y_n * 0 + NI_New1_z_n * 1;

						//������� � ����� ������ ���������
						NI_New1_z_n = NI_New1_z_n + SC[s].onePoint.z;

						if (NI[i][j].onePointZ > NI_New1_z_n)
							NI[i][j].onePointZ = NI_New1_z_n;

						if (NI[i][j].twoPointZ < NI_New1_z_n)
							NI[i][j].twoPointZ = NI_New1_z_n;
					}
				}
			}
		});
    }
	NI.synchronize();
}

Point randomizePoint(double leftLimit, double rightLimit)
{
	Point newPoint;
	
	newPoint.x = leftLimit + (rand() % static_cast<int>(rightLimit - leftLimit + 1));
	newPoint.y = leftLimit + (rand() % static_cast<int>(rightLimit - leftLimit + 1));
	newPoint.z = leftLimit + (rand() % static_cast<int>(rightLimit - leftLimit + 1));

	return newPoint;
}

Capsule generateCapsule(double minLength, double maxLength, double rightLimit)
{
	if (minLength < 0.0)
		minLength = 0.0;

	if (maxLength < 0.0)
		maxLength = 0.0;

	if (minLength > 200000.0)
		minLength = 200000.0;

	if (maxLength > 200000.0)
		maxLength = 200000.0;

	Capsule newCapsule;

	newCapsule.onePoint = randomizePoint(0, rightLimit);
	Point startPoint = newCapsule.onePoint;
	
	while (true)
	{
		newCapsule.twoPoint = randomizePoint(0, rightLimit);

		Point endPoint = newCapsule.twoPoint;

		double length = abs(
							sqrt(
								(startPoint.x - endPoint.x) * (startPoint.x - endPoint.x)
							  + (startPoint.y - endPoint.y) * (startPoint.y - endPoint.y)
							  + (startPoint.z - endPoint.z) * (startPoint.z - endPoint.z)
							)
						);
		if (length <= maxLength && length >= minLength)
			return newCapsule;
	}
}

void Calc_0(vector<TwoPoints>& NeedIntersections, double R )
{

	for (int i = 0; i < NeedIntersections.size(); i++)
	{
		if (NeedIntersections[i].onePointZ == 100000 * R)
		{
			NeedIntersections[i].onePointZ = 0;
		}

		if (NeedIntersections[i].twoPointZ == -100000 * R)
		{
			NeedIntersections[i].twoPointZ = 0;
		}
	}

}

int main()
{
	setlocale(LC_CTYPE, "rus");
	int StepX;
	int	StepY;
	double R;
	int list;
	int NumberOfScope; //TODO: ����������� ����

	accelerator defaultDevice(accelerator::default_accelerator);
	accelerator_view defaultView = defaultDevice.default_view;
	wcout << L" Using device : " << defaultDevice.get_description() << endl << endl; //TODO: �������� �������������� ����������

	cout << "���������� �����" << endl;
	cin >> StepX;

	StepY = StepX;

	const int ArraySize = (StepX + 1) * (StepY + 1);

	cout << "������� ������ �����������" << endl;
	cin >> R;

	cout << "����������� ������" << endl;
	cin >> NumberOfScope;

	vector <Capsule> ScopeCenter(NumberOfScope); //����� ������

	for (int i = 0; i < NumberOfScope; i++)
	{
		ScopeCenter[i] = generateCapsule(R, 5*R, StepX);

	}

	vector<TwoPoints> NeedIntersections(ArraySize);

	for (int i = 0; i < ArraySize; i++)
	{
		NeedIntersections[i].onePointZ = 100000 * R;
		NeedIntersections[i].twoPointZ = -100000 * R;
	}

	cout << "GPU_v3" << std::endl;
	unsigned int start_time = clock(); //��������� �����

	GPUCalc_3(StepX, StepX, R, ScopeCenter, NeedIntersections);
	Calc_0(NeedIntersections, R);

	unsigned int end_time = clock(); // �������� �����
	unsigned int search_time = end_time - start_time; // ������� �����
	cout << search_time << " ms." << endl;


	NeedIntersections.clear();

	system("pause");
}
