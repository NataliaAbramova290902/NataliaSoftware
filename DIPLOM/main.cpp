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
	double x;
	double y;
	double z;
};


struct TwoPoints
{
	Point onePoint;
	Point twoPoint;
};

struct IntersectionPoints
{
	vector<TwoPoints> twoPoints;
};

struct Line
{
		Point onePoint;
		Point twoPoint;
		//vector<IntersectionPoint> IntersectionPoints;
};

void SimpleSqrtGPUCalc(std::vector<Line>& Lines, double R, vector<Point>& ScopeCenter, vector< vector<TwoPoints>> NeedIntersections)
{
	int N = Lines.size();
	array_view<const Line, 1> L(N, Lines);
	array_view<const Point, 1> SC(ScopeCenter.size(), ScopeCenter);

	array_view<TwoPoints, 2> NI(ScopeCenter.size(), N, NeedIntersections);
	NI.discard_data();
	
	parallel_for_each(NI.extent, [=](index<2> idx) restrict(amp)
		//parallel_for(size_t(0), Lines.size(), [&](int idx)
		//for (int idx = 0; idx < Lines.size(); idx++)
	{

        //IntersectionPoints IP;
		
		
		
		
		for (int j = 0; j < N; j++)
		{

			double A = (L(j).twoPoint.x - L(j).onePoint.x) * (L(j).twoPoint.x - L(j).onePoint.x) + (L(j).twoPoint.y - L(j).onePoint.y) * (L(j).twoPoint.y - L(j).onePoint.y) + (L(j).twoPoint.z - L(j).onePoint.z) * (L(j).twoPoint.z - L(j).onePoint.z);
			double B = 2 * ((L(j).twoPoint.x - L(j).onePoint.x) * (L(j).onePoint.x - SC(idx).x) + (L(j).twoPoint.y - L(idx).onePoint.y) * (L(j).onePoint.y - SC(idx).y) + (L(j).twoPoint.z - L(j).onePoint.z) * (L(j).onePoint.z - SC(idx).z));
			double C = (L(j).onePoint.x - SC(idx).x) * (L(j).onePoint.x - SC(idx.x) + (L(j).onePoint.y - SC(idx).y) * (L(j).onePoint.y - SC(idx).y) + (L(j).onePoint.z - SC(idx).z) *(L(j).onePoint.z - SC(idx).z) - (R * R);

			double D = (B * B) - (4 * A * C);

			if (D > 0)
			{
				double SQ = fast_math::sqrt(D);
				double t1 = (-B + SQ) / (2 * A);
				double t2 = (-B - SQ) / (2 * A);

				//Первая точка пересечения
				double x1 = L(j).onePoint.x + t1 * (L(j).twoPoint.x - L(j).onePoint.x);
				double y1 = L(j).onePoint.y + t1 * (L(j).twoPoint.y - L(j).onePoint.y);
				double z1 = L(j).onePoint.z + t1 * (L(j).twoPoint.z - L(j).onePoint.z);

				//Вторая точка пересечения
				double x2 = L(j).onePoint.x + t2 * (L(j).twoPoint.x - L(j).onePoint.x);
				double y2 = L(j).onePoint.y + t2 * (L(j).twoPoint.y - L(j).onePoint.y);
				double z2 = L(j).onePoint.z + t2 * (L(j).twoPoint.z - L(j).onePoint.z);

				Point a1; 
				a1.x = x1;
				a1.y = y1;
				a1.z = z1;
				//IP.Points.push_back(std::move(a1));

				
				NI[idx].onePoint


				Point a2;
				a2.x = x2;
				a2.y = y2;
				a2.z = z2;
				//IP.Points = x1;

				//IP.Points.push_back(std::move(a2));
				//NI[idx] = IP;
				//NI.push_back(IP2);
			}

			else
			{
				if (D == 0)
				{
					double SQ = fast_math::sqrt(D);
					double t1 = (-B + SQ) / (2 * A);
					//Первая точка пересечения
					double x1 = L(j).onePoint.x + t1 * (L(j).twoPoint.x - L(j).onePoint.x);
					double y1 = L(j).onePoint.y + t1 * (L(j).twoPoint.y - L(j).onePoint.y);
					double z1 = L(j).onePoint.z + t1 * (L(j).twoPoint.z - L(j).onePoint.z);
					
					Point a1;
					a1.x = x1;
					a1.y = y1;
					a1.z = z1;
					//IP.Points.data.push_back(a1);
				}
			}
		}
		//NI[idx] = IP;
   });
	NI.synchronize();
}


void CPUCalc(std::vector<Line>& Lines, double R, vector<Point>& ScopeCenter)
{
	vector<IntersectionPoints> NeedIntersections;

	for (int j = 0; j < ScopeCenter.size(); j++) 
	{
		IntersectionPoints IP;
		for (int i = 0; i < Lines.size(); i++) 
		{

			double A = (Lines.at(i).twoPoint.x - Lines.at(i).onePoint.x) * (Lines.at(i).twoPoint.x - Lines.at(i).onePoint.x) + (Lines.at(i).twoPoint.y - Lines.at(i).onePoint.y) * (Lines.at(i).twoPoint.y - Lines.at(i).onePoint.y) + (Lines.at(i).twoPoint.z - Lines.at(i).onePoint.z) * (Lines.at(i).twoPoint.z - Lines.at(i).onePoint.z);
			double B = 2 * ((Lines.at(i).twoPoint.x - Lines.at(i).onePoint.x) * (Lines.at(i).onePoint.x - ScopeCenter.at(j).x) + (Lines.at(i).twoPoint.y - Lines.at(i).onePoint.y) * (Lines.at(i).onePoint.y - ScopeCenter.at(j).y) + (Lines.at(i).twoPoint.z - Lines.at(i).onePoint.z) * (Lines.at(i).onePoint.z - ScopeCenter.at(j).z));
			double C = (Lines.at(i).onePoint.x - ScopeCenter.at(j).x) * (Lines.at(i).onePoint.x - ScopeCenter.at(j).x) + (Lines.at(i).onePoint.y - ScopeCenter.at(j).y) * (Lines.at(i).onePoint.y - ScopeCenter.at(j).y) + (Lines.at(i).onePoint.z - ScopeCenter.at(j).z) *(Lines.at(i).onePoint.z - ScopeCenter.at(j).z) - (R * R);

			double D = (B * B) - (4 * A * C);
			double SQ = fast_math::sqrt(D);
			if (D > 0)
			{
				double t1 = (-B + SQ) / (2 * A);
				double t2 = (-B - SQ) / (2 * A);

				//Первая точка пересечения
				double x1 = Lines.at(i).onePoint.x + t1 * (Lines.at(i).twoPoint.x - Lines.at(i).onePoint.x);
				double y1 = Lines.at(i).onePoint.y + t1 * (Lines.at(i).twoPoint.y - Lines.at(i).onePoint.y);
				double z1 = Lines.at(i).onePoint.z + t1 * (Lines.at(i).twoPoint.z - Lines.at(i).onePoint.z);

				//Вторая точка пересечения
				double x2 = Lines.at(i).onePoint.x + t2 * (Lines.at(i).twoPoint.x - Lines.at(i).onePoint.x);
				double y2 = Lines.at(i).onePoint.y + t2 * (Lines.at(i).twoPoint.y - Lines.at(i).onePoint.y);
				double z2 = Lines.at(i).onePoint.z + t2 * (Lines.at(i).twoPoint.z - Lines.at(i).onePoint.z);


				Point a1;
				a1.x = x1;
				a1.y = y1;
				a1.z = z1;
				IP.Points.push_back(std::move(a1));

				Point a2;
				a2.x = x2;
				a2.y = y2;
				a2.z = z2;
				IP.Points.push_back(std::move(a2));

				NeedIntersections.push_back(IP);
			}

			else
			{
				if (D == 0)
				{
					double t1 = (-B + SQ) / (2 * A);
					//Первая точка пересечения
					double x1 = Lines.at(i).onePoint.x + t1 * (Lines.at(i).twoPoint.x - Lines.at(i).onePoint.x);
					double y1 = Lines.at(i).onePoint.y + t1 * (Lines.at(i).twoPoint.y - Lines.at(i).onePoint.y);
					double z1 = Lines.at(i).onePoint.z + t1 * (Lines.at(i).twoPoint.z - Lines.at(i).onePoint.z);
					IntersectionPoints IP;
					IP.p1.x = x1;
					IP.p1.y = y1;
					IP.p1.z = z1;
					//IP1.Scope = j;
					IP.p2.x = NULL;
					IP.p2.y = NULL;
					IP.p2.z = NULL;

					NeedIntersections.push_back(IP);
				}
				else
				{
					IntersectionPoints IP;
					IP.p1.x = NULL;
					IP.p1.y = NULL;
					IP.p1.z = NULL;

					IP.p2.x = NULL;
					IP.p2.y = NULL;
					IP.p2.z = NULL;

					NeedIntersections.push_back(IP);
				}
			}
		}

		//for (int i = 0; i < NeedIntersections.size(); i++)
		//{
			//Lines[i].IntersectionPoints.push_back(NeedIntersections.at(i));
		//}

		//NeedIntersections.clear();

	}
};



int main()
{
	setlocale(LC_CTYPE, "rus");
	int StepX;
	int	StepY;
	double R;
	int NumberOfScope; //колличество сфер
	//Point ScopeCenter[2]; //Центры сфер
	vector <Point> ScopeCenter;
	accelerator defaultDevice(accelerator::default_accelerator);
	accelerator_view defaultView = defaultDevice.default_view;
	wcout << L" Using device : " << defaultDevice.get_description() << endl << endl; //Название использующейся видеокарты
	
	cout << "Количество шагов по оси X" << endl;
	cin >> StepX;
	
	cout << "Количество шагов по оси Y" << endl;
	cin >> StepY;
	
	const int ArraySize = StepX + StepY + 2;
	
	cout << "Введите радиус" << endl;
	cin >> R;
	
	//ifstream F;
	//F.open("S.txt", ios::in); // окрываем файл для чтения
	//if (F)
	//{
		//F >> NumberOfScope; //колличество сфер
		cout << "Колличество сфер " ;
		cin >> NumberOfScope;
		//cout << "Координаты центров сфер:\n" << endl;
		for (int i = 0; i < NumberOfScope; i++)
		{
			//ScopeCenter(NumberOfScope); //Центры сфер
			Point C;
			
			C.x = 100;
			C.y = 100;
			C.z = 0;

			ScopeCenter.push_back(C);  // добавляем введенные числа
			//cout << ScopeCenter[ScopeCenter.size() - 1].x << "\t" << ScopeCenter[ScopeCenter.size() - 1].y << "\t" << ScopeCenter[ScopeCenter.size() - 1].z << "\t";
			//cout << "\n";

		//}
		//cout << ScopeCenter.size() << endl;
	}
	int NN = ScopeCenter.size();
	vector<IntersectionPoint> NeedIntersections(NN);
	//Координаты угла левого нижнего
	Point U;
	Point U2;
	U.x = 0.0;
	U.y = 0.0;
	U.z = 0.0;
	//Координаты угла правого верхнего угла
	U2.x = 300.0;
	U2.y = 300.0;
	U2.z = 0.0;

	//кординаты линий

	vector<Line> Lines(ArraySize); 		
	for (int i = 0; i < StepX+1; i++)
	{
        Lines[i].onePoint.x = ((U2.x - U.x) / StepX) * i; //координаты X первой точки линии параллельных оси У
		Lines[i].onePoint.y = U.y; //координаты Y первой точки линии параллельных оси У
		Lines[i].onePoint.z = U.z; //координаты Z первой точки линии параллельных оси У
		Lines[i].twoPoint.x = ((U2.x - U.x) / StepX) * i;//координаты X второй точки линии параллельных оси У
		Lines[i].twoPoint.y = U2.y; //координаты Y второй точки линии параллельных оси У
		Lines[i].twoPoint.z = U.z; //координаты Z второй точки линии параллельных оси У
	}

	for (int i = StepX + 1; i < ArraySize; i++)
	{
		Lines[i].onePoint.y = ((U2.y - U.y) / StepY) * (i- StepX - 1); //координаты У первой точки линии параллельных оси Х
		Lines[i].onePoint.x = U.x; //координаты Х первой точки линии параллельных оси Х
		Lines[i].onePoint.z = U.z; //координаты Z первой точки линии параллельных оси X
		Lines[i].twoPoint.y = ((U2.y - U.y) / StepY) * (i - StepX - 1);//координаты У второй точки линии параллельных оси Х
		Lines[i].twoPoint.x = U2.x;//координаты Х второй точки линии параллельных оси Х
		Lines[i].twoPoint.z = U2.z; //координаты Z второй точки линии параллельных оси X
	}

	/*cout << "Точки линий" << endl;
	cout << "\n";
	for (int i = 0; i < ArraySize; i++)
	{
		cout << "Линия " << i+1 << endl;
		cout << Lines[i].onePoint.x << "\t" << Lines[i].onePoint.y << "\t" << Lines[i].onePoint.z << "\t" << endl;
		cout << Lines[i].twoPoint.x << "\t" << Lines[i].twoPoint.y << "\t" << Lines[i].twoPoint.z << "\t" << endl;
		cout << "\n";

	}*/

	cout << "GPU" << std::endl << std::endl;
	unsigned int start_time = clock(); //начальное время

	SimpleSqrtGPUCalc(Lines, R, ScopeCenter, NeedIntersections);

	unsigned int end_time = clock(); // конечное время
	unsigned int search_time = end_time - start_time; // искомое время
	cout << search_time << " ms." << endl;

	cout << "СPU" << endl;
	unsigned int start_time2 = clock(); //начальное время

	CPUCalc(Lines, R, ScopeCenter);

	unsigned int end_time2 = clock(); // конечное время
	unsigned int search_time2 = end_time2 - start_time2; // искомое время
	cout << search_time2 << " ms." << endl;

	/*for (int i = 0; i < Lines.size(); i++)
	{		
		for (int j = 0; j < Lines.at(i).IntersectionPoints.size(); j++)
		{	
			cout << "Линия " << i + 1 << ", Сфера " << Lines[i].IntersectionPoints[j].Scope + 1 << ", X = " << Lines[i].IntersectionPoints[j].x << ", Y = " << Lines.at(i).IntersectionPoints.at(j).y << ", Z = " << Lines.at(i).IntersectionPoints.at(j).z << "\n";
		}
	}*/
		
	system("pause");
}