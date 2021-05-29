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

//Структура для хранения точек капсул
struct Capsule
{
	static double radius;
	Point onePoint;
	Point twoPoint;
	
};


void GPUCalc_3(int StepX, int StepY, double R, vector<Capsule>& ScopeCenter, vector<TwoPoints>& NeedIntersections)
{
	
	array_view<const Point, 1 > SC(ScopeCenter.size(), ScopeCenter);
	array_view<TwoPoints, 2> NI((StepX + 1) , (StepY + 1), NeedIntersections);
	NI.discard_data();
	
	for (int s = 0; s < ScopeCenter.size(); s++)
	//parallel_for_each(SC.extent, [=](index<1> idx) restrict(amp)
	{
		int NI2y = (fast_math::ceil(SC[s].y + R) + 1);
		int NI1y = (fast_math::floor(SC[s].y - R));

		int NI2x = (fast_math::ceil(SC[s].x + R) + 1);
		int NI1x = (fast_math::floor(SC[s].x - R));

		concurrency::extent<2> e(NI2y - NI1y, NI2x - NI1x);
		parallel_for_each(e, [=](index<2> idx) restrict(amp)
			{
				int j = (NI1x + idx[1]);
				int i = (NI1y + idx[0]);
				double A = ((-100000 * R) - (100000 * R)) * ((-100000 * R) - (100000 * R));
				double B = 2 * (((-100000 * R) - (100000 * R)) * ((100000 * R) - SC[s].z));
				double C = ((j - SC[s].x) * (j - SC[s].x) + (i - SC[s].y) * (i - SC[s].y) + ((100000 * R) - SC[s].z) * ((100000 * R) - SC[s].z)) - (R * R);

				double D = (B * B) - 4 * A * C;

				if (D > 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					double t2 = (-B - SQ) / (2 * A);

					//Первая точка пересечения
					double x1 = j;
					double y1 = i;
					double z1 = (100000 * R) + t1 * ((-100000 * R) - (100000 * R));

					//Вторая точка пересечения
					double x2 = j;
					double y2 = i;
					double z2 = (100000 * R) + t2 * ((-100000 * R) - (100000 * R));

					if (z1 < z2)
					{
						if(NI[i][j].onePointZ > z1)
							NI[i][j].onePointZ = z1;

						if(NI[i][j].twoPointZ < z2)
							NI[i][j].twoPointZ = z2;
					}
					else
					{
						if (NI[i][j].onePointZ > z2)
							NI[i][j].onePointZ = z2;

						if (NI[i][j].twoPointZ < z1)
							NI[i][j].twoPointZ = z1;
					}
				}

				else if (D == 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					//Первая точка пересечения
					double x1 = j;
					double y1 = i;
					double z1 = (100000 * R) + t1 * ((-100000 * R) - (100000 * R));

					double x2 = j;
					double y2 = i;
					double z2 = (100000 * R) + t1 * ((-100000 * R) - (100000 * R));

					if (NI[i][j].onePointZ > z1)
						NI[i][j].onePointZ = z1;

					if (NI[i][j].twoPointZ < z1)
						NI[i][j].twoPointZ = z1;

				}
			//}

		});
    }
	NI.synchronize();
	
	
}


void CPUCalc(double R, vector<Capsule>& ScopeCenter, vector<TwoPoints>& NeedIntersections, int StepX)
{
	for (int k = 0; k < ScopeCenter.size(); k++)
	{ 
		for (int i = 0; i < StepX; i++)
		{
			for (int j = 0; j < StepX; j++)
			{
				double A = ((-100000 * R) - (100000 * R)) * ((-100000 * R) - (100000 * R));
				double B = 2 * (((-100000 * R) - (100000 * R)) * ((100000 * R) - ScopeCenter[k].z));
				double C = ((j - ScopeCenter[k].x) * (j - ScopeCenter[k].x) + (i - ScopeCenter[k].y) * (i - ScopeCenter[k].y) + ((100000 * R) - ScopeCenter[k].z) * ((100000 * R) - ScopeCenter[k].z)) - (R * R);


				double D = (B * B) - 4 * A * C;

				if (D > 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					double t2 = (-B - SQ) / (2 * A);

					//TODO:Первая точка пересечения
					double z1 = (100000 * R) + t1 * ((-100000 * R) - (100000 * R));

					//TODO:Вторая точка пересечения
					double z2 = (100000 * R) + t2 * ((-100000 * R) - (100000 * R));

					if (z1 < z2)
					{
						double NPPP = NeedIntersections[i * (StepX + 1) + j].onePointZ;
						if(NeedIntersections[i * (StepX + 1) + j].onePointZ > z1)
							NeedIntersections[i * (StepX + 1) + j].onePointZ = z1;

						if(NeedIntersections[i * (StepX + 1) + j].twoPointZ < z2)
							NeedIntersections[i * (StepX + 1) + j].twoPointZ = z2;
					}
					else
					{
						if (NeedIntersections[i * (StepX + 1) + j].onePointZ > z2)
							NeedIntersections[i * (StepX + 1) + j].onePointZ = z2;

						if (NeedIntersections[i * (StepX + 1) + j].twoPointZ < z1)
							NeedIntersections[i * (StepX + 1) + j].twoPointZ = z1;
					}
				}

				else if (D == 0)
				{
					double SQ = sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					//TODO: Первая точка пересечения
					double z1 = (100000 * R) + t1 * ((-100000 * R) - (100000 * R));

					if (NeedIntersections[i * (StepX + 1) + j].onePointZ > z1)
						NeedIntersections[i * (StepX + 1) + j].onePointZ = z1;

					if (NeedIntersections[i * (StepX + 1) + j].twoPointZ < z1)
						NeedIntersections[i * (StepX + 1) + j].twoPointZ = z1;


				}
			}

		}
	}
};

/*void random(Capsule Center, double R, int StepX, int StepY)
{
	if((Center.twoPoint.x - Center.onePoint.x) * (Center.twoPoint.x - Center.onePoint.x) + (Center.twoPoint.y - Center.onePoint.y) * (Center.twoPoint.y - Center.onePoint.y) + (Center.twoPoint.z - Center.onePoint.z) * (Center.twoPoint.z - Center.onePoint.z) - (R * R) < 0)
	{ 
		
		Center.twoPoint.x = rand() % StepX;
		Center.twoPoint.y = rand() % StepY;
		Center.twoPoint.z = rand();
		random(Center, R, StepX, StepY);
	}
	return  Center;
}*/


Point randomizePoint(double leftLimit, double rightLimit)
{
	Point newPoint;
	
	newPoint.x = leftLimit + (rand() % static_cast<int>(rightLimit - leftLimit + 1));
	newPoint.y = leftLimit + (rand() % static_cast<int>(rightLimit - leftLimit + 1));
	newPoint.z = leftLimit + (rand() % static_cast<int>(rightLimit - leftLimit + 1));

	return newPoint;
}

Capsule generateCapsule(double minLength, double maxLength)
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

	newCapsule.onePoint = randomizePoint(-100000.0, 100000.0);
	Point startPoint = newCapsule.onePoint;
	
	while (true)
	{
		newCapsule.twoPoint = randomizePoint(-100000.0, 100000.0);

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


int main()
{
	Capsule::radius = 20;
	
	setlocale(LC_CTYPE, "rus");
	int StepX;
	int	StepY;
	double R;
	int list;
	int NumberOfScope; //TODO: колличество сфер

	accelerator defaultDevice(accelerator::default_accelerator);
	accelerator_view defaultView = defaultDevice.default_view;
	wcout << L" Using device : " << defaultDevice.get_description() << endl << endl; //TODO: Название использующейся видеокарты

	cout << "Количество шагов по оси X" << endl;
	cin >> StepX;

	cout << "Количество шагов по оси Y" << endl;
	cin >> StepY;

	const int ArraySize = (StepX + 1) * (StepY + 1);

	cout << "Введите радиус инструмента" << endl;
	cin >> R;
	int max_c = (floor(StepX / (2*(R + 0.1))) * (floor(StepY / (2*(R + 0.1)))));

	cout << "Колличество капсул" /*(максимальное колиество " << max_c << ")"*/ << endl;
	cin >> NumberOfScope;

    

	vector <Capsule> ScopeCenter(NumberOfScope); //Точки капсул

	for (int i = 0; i < NumberOfScope; i++)
	{
		ScopeCenter[i].onePoint.x = rand() % StepX;
		ScopeCenter[i].onePoint.y = rand() % StepY;
		ScopeCenter[i].onePoint.z = rand();

		ScopeCenter[i].twoPoint.x = rand() % StepX;
		ScopeCenter[i].twoPoint.y = rand() % StepY;
		ScopeCenter[i].twoPoint.z = rand();

		//random(ScopeCenter[i]);
	}

	vector<TwoPoints> NeedIntersections(ArraySize);


	cout << "GPU_v3" << std::endl;
	unsigned int start_time = clock(); //начальное время

	GPUCalc_3(StepX, StepX, R, ScopeCenter, NeedIntersections);

	unsigned int end_time = clock(); // конечное время
	unsigned int search_time = end_time - start_time; // искомое время
	cout << search_time << " ms." << endl;

	NeedIntersections.clear();

	vector<TwoPoints>().swap(NeedIntersections);

	vector<TwoPoints> NeedIntersections4;

	NeedIntersections4.resize(ArraySize);


	cout << "СPU" << endl;
	unsigned int start_time4 = clock(); //начальное время

	CPUCalc(R, ScopeCenter, NeedIntersections4, StepX);

	unsigned int end_time4 = clock(); // конечное время
	unsigned int search_time4 = end_time4 - start_time4; // искомое время
	cout << search_time4 << " ms." << endl;

	system("pause");
}
