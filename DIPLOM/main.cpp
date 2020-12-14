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
	double x = 0;
	double y = 0;
	double z = 0;
	//bool valide = true;   //TODO: Флаг валидности точки, типа в теории мы все равно должны создать экземпляр - прототип точки, чтоб отроботал расчет
};

//TODO: Структура для храенения точек пересечения, типа их не может быть больше двух
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

void GPUCalc(int StepX, int StepY, double R, vector<Point>& ScopeCenter, vector<TwoPoints>& NeedIntersections)
{
	//int N = NeedIntersections.size();
	//array_view<const Line, 1 > L((StepX+1) * (StepY+1), Lines);
	array_view<const Point, 1 > SC(ScopeCenter.size(), ScopeCenter);
	array_view<TwoPoints, 1> NI((StepX + 1) * (StepY + 1), NeedIntersections);
	NI.discard_data();
	
	
	parallel_for_each(SC.extent, [=](index<1> idx) restrict(amp)
	{
		for (int i = fast_math::floor(SC[idx].y - R); i < fast_math::floor(SC[idx].y + R) + 1; i++)
		{
			for (int j = fast_math::floor(SC[idx].x - R); j < fast_math::floor(SC[idx].x + R) + 1; j++)
			{
				double A = ((-R - 0.1) - (R + 0.1)) * ((-R - 0.1) - (R + 0.1));
				double B = 2 * (((-R - 0.1) - (R + 0.1)) * ((R + 0.1) - SC[idx].z));
				double C = ((j - SC[idx].x) * (j - SC[idx].x) + (i - SC[idx].y) * (i - SC[idx].y) + ((R + 0.1) - SC[idx].z) * ((R + 0.1) - SC[idx].z)) - (R * R);

				double D = (B * B) - 4 * A * C;

				if (D > 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					double t2 = (-B - SQ) / (2 * A);

					//Первая точка пересечения
					double x1 = j;
					double y1 = i;
					double z1 = (R + 0.1) + t1 * ((-R - 0.1) - (R + 0.1));

					//Вторая точка пересечения
					double x2 = j;
					double y2 = i;
					double z2 = (R + 0.1) + t2 * ((-R - 0.1) - (R + 0.1));

					///NI[i * (StepX + 1) + j].onePoint.x = x1;
					//NI[i * (StepX + 1) + j].onePoint.y = y1;
					NI[i * (StepX + 1) + j].onePointZ = z1;
					//NI[i * (StepX + 1) + j].onePoint.valide = true;

					//NI[i * (StepX + 1) + j].twoPoint.x = x2;
					//NI[i * (StepX + 1) + j].twoPoint.y = y2;
					NI[i * (StepX + 1) + j].twoPointZ = z2;
					//NI[i * (StepX + 1) + j].twoPoint.valide = true;
				}

				else if (D == 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					//Первая точка пересечения
					double x1 = j;
					double y1 = i;
					double z1 = (R + 0.1) + t1 * ((-R - 0.1) - (R + 0.1));

					double x2 = j;
					double y2 = i;
					double z2 = (R + 0.1) + t1 * ((-R - 0.1) - (R + 0.1));

					//NI[i * (StepX + 1) + j].onePoint.x = x1;
					//NI[i * (StepX + 1) + j].onePoint.y = y1;
					NI[i * (StepX + 1) + j].onePointZ = z1;
					//NI[i * (StepX + 1) + j].onePoint.valide = true;

					NI[i * (StepX + 1) + j].twoPointZ = NI[i * (StepX + 1) + j].onePointZ;

				}
			}

		}
    });
	NI.synchronize();
	
	
}


void CPUCalc(double R, vector<Point>& ScopeCenter, vector<TwoPoints>& NeedIntersections, int StepX)
{
	for (int k = 0; k < ScopeCenter.size(); k++)
	{ 
		for (int i = floor(ScopeCenter[k].y - R); i < floor(ScopeCenter[k].y + R)+1; i++)
		{
			for (int j = floor(ScopeCenter[k].x - R); j < floor(ScopeCenter[k].x + R)+1; j++)
			{
				double A = ((-R - 0.1) - (R + 0.1)) * ((-R - 0.1) - (R + 0.1));
				double B = 2 * (((-R - 0.1) - (R + 0.1)) * ((R + 0.1) - ScopeCenter[k].z));
				double C = ((j - ScopeCenter[k].x) * (j - ScopeCenter[k].x) + (i - ScopeCenter[k].y) * (i - ScopeCenter[k].y) + ((R + 0.1) - ScopeCenter[k].z) * ((R + 0.1) - ScopeCenter[k].z)) - (R * R);


				double D = (B * B) - 4 * A * C;

				if (D > 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					double t2 = (-B - SQ) / (2 * A);

					//TODO:Первая точка пересечения
					double x1 = j;
					double y1 = i;
					double z1 = (R + 0.1) + t1 * ((-R - 0.1) - (R + 0.1));

					//TODO:Вторая точка пересечения
					double x2 = j;
					double y2 = i;
					double z2 = (R + 0.1) + t2 * ((-R - 0.1) - (R + 0.1));


					double a1;
					//a1.x = x1;
					//a1.y = y1;
					a1 = z1;

					NeedIntersections[i * (StepX + 1) + j].onePointZ = a1;

					double a2;
					//a2.x = x2;
					//a2.y = y2;
					a2 = z2;
					NeedIntersections[i * (StepX + 1) + j].twoPointZ = a2;
				}

				else if (D == 0)
				{
					double SQ = sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					//TODO: Первая точка пересечения
					double x1 = j;
					double y1 = i;
					double z1 = (R + 0.1) + t1 * ((-R - 0.1) - (R + 0.1));

					double a1;
					//a1.x = x1;
					//a1.y = y1;
					a1 = z1;
					//a1.valide = true;

					NeedIntersections[i * (StepX + 1) + j].onePointZ = a1;

					NeedIntersections[i * (StepX + 1) + j].twoPointZ = NeedIntersections[i * (StepX + 1) + j].onePointZ;

				}
			}

		}
	}
};



int main()
{
	setlocale(LC_CTYPE, "rus");
	int StepX;
	int	StepY;
	double R;
	int NumberOfScope; //TODO: колличество сфер

	accelerator defaultDevice(accelerator::default_accelerator);
	accelerator_view defaultView = defaultDevice.default_view;
	wcout << L" Using device : " << defaultDevice.get_description() << endl << endl; //TODO: Название использующейся видеокарты

	cout << "Количество шагов по оси X" << endl;
	cin >> StepX;

	cout << "Количество шагов по оси Y" << endl;
	cin >> StepY;

	const int ArraySize = (StepX + 1) * (StepY + 1);

	cout << "Введите радиус" << endl;
	cin >> R;
	int max_c = (floor(StepX / (2*(R + 0.1))) * (floor(StepY / (2*(R + 0.1)))));

	cout << "Колличество сфер (максимальное колиество " << max_c << ")" << endl;
	cin >> NumberOfScope;


	int stringCount= floor(NumberOfScope / (floor(StepY / 2 / (R + 0.1))));
	int rowCount = floor(StepX / 2 / (R + 0.1));

	vector <Point> ScopeCenter(NumberOfScope); //Центры сфер	
	for (int i = 0; i < stringCount; ++i)
	{
		for (int j = 0; j < rowCount; j++)
		{
			ScopeCenter[i * floor(StepY / 2 /(R + 0.1)) + j].x = (R + 0.1) * (2 * j + 1);
			ScopeCenter[i * floor(StepY / 2 /(R + 0.1)) + j].y = (R + 0.1) * (2 * i + 1);
			ScopeCenter[i * floor(StepY / 2 /(R + 0.1)) + j].z = 0;
		}
	}
	int n = floor(StepY / 2 / (R + 0.1));
	if (NumberOfScope % n != 0)
	{
		for (int i = (floor(NumberOfScope / floor(StepY / 2 / (R + 0.1)))) * floor(StepX / 2 / (R + 0.1)); i < (floor(NumberOfScope / floor(StepY / 2 / (R + 0.1)))) * floor(StepX / 2 / (R + 0.1)) + 1; i++)
		{
			for (int j = 0; j < (NumberOfScope - (floor(NumberOfScope / floor(StepY / 2 / (R + 0.1)))) * floor(StepX / 2 / (R + 0.1))); j++)
			{
				ScopeCenter[i + j].x = (R + 0.1) * (2 * j + 1);
				ScopeCenter[i + j].y = (R + 0.1) * (floor(NumberOfScope / (floor(StepY / 2 / (R + 0.1)-2))));
				ScopeCenter[i + j].z = 0.0;
			}
		}
	}
		


	int NN = ScopeCenter.size();
	vector<TwoPoints> NeedIntersections(ArraySize);

	cout << "GPU" << std::endl;
	unsigned int start_time = clock(); //начальное время

	GPUCalc(StepX, StepX, R, ScopeCenter, NeedIntersections);

	unsigned int end_time = clock(); // конечное время
	unsigned int search_time = end_time - start_time; // искомое время
	cout << search_time << " ms." << endl;



	NeedIntersections.clear();

	NeedIntersections = vector<TwoPoints>(ArraySize);


	cout << "СPU" << endl;
	unsigned int start_time2 = clock(); //начальное время

	CPUCalc(R, ScopeCenter, NeedIntersections, StepX);

	unsigned int end_time2 = clock(); // конечное время
	unsigned int search_time2 = end_time2 - start_time2; // искомое время
	cout << search_time2 << " ms." << endl;
	system("pause");
}
