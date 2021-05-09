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


void GPUCalc_1(int StepX, int StepY, double R, vector<Point>& ScopeCenter, vector<TwoPoints>& NeedIntersections)
{
	array_view<const Point, 1 > SC(ScopeCenter.size(), ScopeCenter);
	array_view<TwoPoints, 2> NI((StepX + 1), (StepY + 1), NeedIntersections);
	NI.discard_data();


	parallel_for_each(SC.extent, [=](index<1> idx) restrict(amp)
	{
		for (int i = 0; i < StepX + 1; i++)
		{
			for (int j = 0; j < StepY + 1; j++)
			{
				double A = ((-100000 * R) - (100000 * R)) * ((-100000 * R) - (100000 * R));
				double B = 2 * (((-100000 * R) - (100000 * R)) * ((100000 * R) - SC[idx].z));
				double C = ((j - SC[idx].x) * (j - SC[idx].x) + (i - SC[idx].y) * (i - SC[idx].y) + ((100000 * R) - SC[idx].z) * ((100000 * R) - SC[idx].z)) - (R * R);

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
						if (NI[i][j].onePointZ > z1)
							NI[i][j].onePointZ = z1;

						if (NI[i][j].twoPointZ < z2)
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
			}

		}
	});
	NI.synchronize();


}

void GPUCalc_2(int StepX, int StepY, double R, vector<Point>& ScopeCenter, vector<TwoPoints>& NeedIntersections, int NumberOfScope)
{
	array_view<const Point, 1 > SC(ScopeCenter.size(), ScopeCenter);
	array_view<TwoPoints, 2> NI((StepX + 1), (StepY + 1), NeedIntersections);
	NI.discard_data();


	parallel_for_each(NI.extent, [=](index<2> idx) restrict(amp)
	{
			for (int c = 0; c < NumberOfScope; c++)
			{
				double A = ((-100000 * R) - (100000 * R)) * ((-100000 * R) - (100000 * R));
				double B = 2 * (((-100000 * R) - (100000 * R)) * ((100000 * R) - SC[c].z));
				double C = ((idx[1] - SC[c].x) * (idx[1] - SC[c].x) + (idx[0] - SC[c].y) * (idx[0] - SC[c].y) + ((100000 * R) - SC[c].z) * ((100000 * R) - SC[c].z)) - (R * R);

				double D = (B * B) - 4 * A * C;

				if (D > 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					double t2 = (-B - SQ) / (2 * A);

					//Первая точка пересечения
					double x1 = idx[1];
					double y1 = idx[0];
					double z1 = (100000 * R) + t1 * ((-100000 * R) - (100000 * R));

					//Вторая точка пересечения
					double x2 = idx[1];
					double y2 = idx[0];
					double z2 = (100000 * R) + t2 * ((-100000 * R) - (100000 * R));

					if (z1 < z2)
					{
						if (NI[idx].onePointZ > z1)
							NI[idx].onePointZ = z1;

						if (NI[idx].twoPointZ < z2)
							NI[idx].twoPointZ = z2;
					}
					else
					{
						if (NI[idx].onePointZ > z2)
							NI[idx].onePointZ = z2;

						if (NI[idx].twoPointZ < z1)
							NI[idx].twoPointZ = z1;
					}
				}

				else if (D == 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					//Первая точка пересечения
					double x1 = idx[1];
					double y1 = idx[0];
					double z1 = (100000 * R) + t1 * ((-100000 * R) - (100000 * R));

					double x2 = idx[1];
					double y2 = idx[0];
					double z2 = (100000 * R) + t1 * ((-100000 * R) - (100000 * R));

					if (NI[idx].onePointZ > z1)
						NI[idx].onePointZ = z1;

					if (NI[idx].twoPointZ < z1)
						NI[idx].twoPointZ = z1;

				}
			}

	});
	NI.synchronize();


}

void GPUCalc_3(int StepX, int StepY, double R, vector<Point>& ScopeCenter, vector<TwoPoints>& NeedIntersections)
{
	array_view<const Point, 1 > SC(ScopeCenter.size(), ScopeCenter);
	array_view<TwoPoints, 2> NI((StepX + 1) , (StepY + 1), NeedIntersections);
	NI.discard_data();
	
	
	parallel_for_each(SC.extent, [=](index<1> idx) restrict(amp)
	{
		for (int i = fast_math::floor(SC[idx].y - R); i < fast_math::ceil(SC[idx].y + R) + 1; i++)
		{
			for (int j = fast_math::floor(SC[idx].x - R); j < fast_math::ceil(SC[idx].x + R) + 1; j++)
			{
				double A = ((-100000 * R) - (100000 * R)) * ((-100000 * R) - (100000 * R));
				double B = 2 * (((-100000 * R) - (100000 * R)) * ((100000 * R) - SC[idx].z));
				double C = ((j - SC[idx].x) * (j - SC[idx].x) + (i - SC[idx].y) * (i - SC[idx].y) + ((100000 * R) - SC[idx].z) * ((100000 * R) - SC[idx].z)) - (R * R);

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
			}

		}
    });
	NI.synchronize();
	
	
}


void CPUCalc(double R, vector<Point>& ScopeCenter, vector<TwoPoints>& NeedIntersections, int StepX)
{
	for (int k = 0; k < ScopeCenter.size(); k++)
	{ 
		for (int i = floor(ScopeCenter[k].y - R); i < ceil(ScopeCenter[k].y + R) + 1; i++)
		{
			for (int j = floor(ScopeCenter[k].x - R); j < ceil(ScopeCenter[k].x + R) + 1; j++)
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



int main()
{
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

	cout << "Введите радиус" << endl;
	cin >> R;
	int max_c = (floor(StepX / (2*(R + 0.1))) * (floor(StepY / (2*(R + 0.1)))));

	cout << "Колличество сфер (максимальное колиество " << max_c << ")" << endl;
	cin >> NumberOfScope;

    

	vector <Point> ScopeCenter(NumberOfScope); //Центры сфер	
	if (NumberOfScope > max_c)
	{
		double list2 = static_cast<double>(NumberOfScope) / static_cast<double>(max_c);
		list = ceil(static_cast<double>(NumberOfScope) / static_cast<double>(max_c));
		int rowCountList = floor(static_cast<double>(StepX) / (2 * (R + 0.1)));
		int stringCountList = floor(static_cast<double>(StepY) / (2 * (R + 0.1)));
		for (int k = 0; k < list - 1; k++)
		{
			for (int i = 0; i < stringCountList; ++i)
			{
				for (int j = 0; j < rowCountList; j++)
				{
					ScopeCenter[(k * rowCountList * stringCountList) + (i * rowCountList + j)].x = (R + 0.1) * (2 * j + 1);
					ScopeCenter[(k * rowCountList * stringCountList) + (i * rowCountList + j)].y = (R + 0.1) * (2 * i + 1);
					ScopeCenter[(k * rowCountList * stringCountList) + (i * rowCountList + j)].z = (R + 0.1) * (2 * k);
				}
			}
		}

		int rowCount = floor(static_cast<double>(StepX) / (2 * (R + 0.1)));
		int stringCount = ceil(static_cast<double>(static_cast<double>(NumberOfScope) - ((static_cast<double>(list) - 1) * static_cast<double>(rowCountList) * static_cast<double>(stringCountList))) / static_cast<double>(rowCount));
		for (int k = list - 1; k < list; k++)
		{
			for (int i = 0; i < stringCount - 1; ++i)
			{
				for (int j = 0; j < rowCount; j++)
				{
					ScopeCenter[(k * rowCountList * stringCountList) + (i * rowCount + j)].x = (R + 0.1) * (2 * j + 1);
					ScopeCenter[(k * rowCountList * stringCountList) + (i * rowCount + j)].y = (R + 0.1) * (2 * i + 1);
					ScopeCenter[(k * rowCountList * stringCountList) + (i * rowCount + j)].z = (R + 0.1) * (2 * k);
				}
			}

			for (int i = stringCount - 1; i < stringCount; i++)
			{
				for (int j = 0; j < (static_cast<double>(NumberOfScope) - ((static_cast<double>(list) - 1) * static_cast<double>(rowCountList) * static_cast<double>(stringCountList))-(static_cast<double>(stringCount) - 1) * static_cast<double>(rowCount)); j++)
				{
					ScopeCenter[(k * rowCountList * stringCountList) + (i * rowCount + j)].x = (R + 0.1) * (2 * j + 1);
					ScopeCenter[(k * rowCountList * stringCountList) + (i * rowCount + j)].y = (R + 0.1) * (2 * i + 1);
					ScopeCenter[(k * rowCountList * stringCountList) + (i * rowCount + j)].z = (R + 0.1) * (2 * k);
				}
			}
		}

	}
	else
	{
		int rowCount = floor(StepX / (2 * (R + 0.1)));
		int stringCount = ceil(static_cast<double>(NumberOfScope) / static_cast<double>(rowCount));
		for (int i = 0; i < stringCount - 1; ++i)
		{
		for (int j = 0; j < rowCount; j++)
			{
				ScopeCenter[i * rowCount + j].x = (R + 0.1) * (2 * j + 1);
				ScopeCenter[i * rowCount + j].y = (R + 0.1) * (2 * i + 1);
				ScopeCenter[i * rowCount + j].z = 0;
			}
		}

		for (int i = stringCount - 1; i < stringCount; i++)
		{
			for (int j = 0; j < (NumberOfScope - ((stringCount - 1) * rowCount)); j++)
			{
				ScopeCenter[i * rowCount + j].x = (R + 0.1) * (2 * j + 1);
				ScopeCenter[i * rowCount + j].y = (R + 0.1) * (2 * i + 1);
				ScopeCenter[i * rowCount + j].z = 0.0;
			}
		}
	}
		
    vector<TwoPoints> N1;
    vector<TwoPoints> N2;

	int NN = ScopeCenter.size();
	vector<TwoPoints> NeedIntersections(ArraySize);

	cout << "GPU_v1" << std::endl;
	unsigned int start_time = clock(); //начальное время

	GPUCalc_1(StepX, StepX, R, ScopeCenter, NeedIntersections);

	unsigned int end_time = clock(); // конечное время
	unsigned int search_time = end_time - start_time; // искомое время
	cout << search_time << " ms." << endl;

	//N1 = NeedIntersections;

	NeedIntersections.clear();

	vector<TwoPoints>().swap(NeedIntersections);

	vector<TwoPoints> NeedIntersections2;

	NeedIntersections2.resize(ArraySize);

	cout << "GPU_v2" << std::endl;
	unsigned int start_time2 = clock(); //начальное время

	GPUCalc_2(StepX, StepX, R, ScopeCenter, NeedIntersections2, NumberOfScope);

	unsigned int end_time2 = clock(); // конечное время
	unsigned int search_time2 = end_time2 - start_time2; // искомое время
	cout << search_time2 << " ms." << endl;


	NeedIntersections.clear();

	vector<TwoPoints>().swap(NeedIntersections);

	vector<TwoPoints> NeedIntersections3;

	NeedIntersections3.resize(ArraySize);

	cout << "GPU_v3" << std::endl;
	unsigned int start_time3 = clock(); //начальное время

	GPUCalc_3(StepX, StepX, R, ScopeCenter, NeedIntersections3);

	unsigned int end_time3 = clock(); // конечное время
	unsigned int search_time3 = end_time3 - start_time3; // искомое время
	cout << search_time3 << " ms." << endl;

	//N1 = NeedIntersections;
	

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

	//N2 = NeedIntersections2;



	system("pause");
}
