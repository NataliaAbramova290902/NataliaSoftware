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
	bool valide = true;   //TODO: Флаг валидности точки, типа в теории мы все равно должны создать экземпляр - прототип точки, чтоб отроботал расчет
};

//TODO: Структура для храенения точек пересечения, типа их не может быть больше двух
struct TwoPoints
{
	Point onePoint;
	Point twoPoint;	
};


struct Line
{
		Point onePoint;
		Point twoPoint;
};

void SimpleSqrtGPUCalc(std::vector<Line>& Lines, double R, vector<Point>& ScopeCenter, vector<TwoPoints>& NeedIntersections)
{
	int N = Lines.size();
	array_view<const Line, 1> L(N, Lines);
	array_view<const Point, 1> SC(ScopeCenter.size(), ScopeCenter);

	array_view<TwoPoints, 2> NI(ScopeCenter.size(), N, NeedIntersections);
	NI.discard_data();
	

	parallel_for_each(NI.extent, [=](index<2> idx) restrict(amp)
	{

			double A = (L(idx[1]).twoPoint.x - L(idx[1]).onePoint.x) * (L(idx[1]).twoPoint.x - L(idx[1]).onePoint.x) + (L(idx[1]).twoPoint.y - L(idx[1]).onePoint.y) * (L(idx[1]).twoPoint.y - L(idx[1]).onePoint.y) + (L(idx[1]).twoPoint.z - L(idx[1]).onePoint.z) * (L(idx[1]).twoPoint.z - L(idx[1]).onePoint.z);
			double B = 2 * ((L(idx[1]).twoPoint.x - L(idx[1]).onePoint.x) * (L(idx[1]).onePoint.x - SC(idx[0]).x) + (L(idx[1]).twoPoint.y - L(idx[0]).onePoint.y) * (L(idx[1]).onePoint.y - SC(idx[0]).y) + (L(idx[1]).twoPoint.z - L(idx[1]).onePoint.z) * (L(idx[1]).onePoint.z - SC(idx[0]).z));
			double C = (L(idx[1]).onePoint.x - SC(idx[0]).x) * (L(idx[1]).onePoint.x - SC(idx[0]).x) + (L(idx[1]).onePoint.y - SC(idx[0]).y) * (L(idx[1]).onePoint.y - SC(idx[0]).y) + (L(idx[1]).onePoint.z - SC(idx[0]).z) *(L(idx[1]).onePoint.z - SC(idx[0]).z) - (R * R);

			double D = (B * B) - (4 * A * C);

			if (D > 0)
			{
				double SQ = fast_math::sqrt(D);
				double t1 = (-B + SQ) / (2 * A);
				double t2 = (-B - SQ) / (2 * A);

				//Первая точка пересечения
				double x1 = L(idx[1]).onePoint.x + t1 * (L(idx[1]).twoPoint.x - L(idx[1]).onePoint.x);
				double y1 = L(idx[1]).onePoint.y + t1 * (L(idx[1]).twoPoint.y - L(idx[1]).onePoint.y);
				double z1 = L(idx[1]).onePoint.z + t1 * (L(idx[1]).twoPoint.z - L(idx[1]).onePoint.z);

				//Вторая точка пересечения
				double x2 = L(idx[1]).onePoint.x + t2 * (L(idx[1]).twoPoint.x - L(idx[1]).onePoint.x);
				double y2 = L(idx[1]).onePoint.y + t2 * (L(idx[1]).twoPoint.y - L(idx[1]).onePoint.y);
				double z2 = L(idx[1]).onePoint.z + t2 * (L(idx[1]).twoPoint.z - L(idx[1]).onePoint.z);

				NI[idx].onePoint.x = x1;
				NI[idx].onePoint.y = y1;
				NI[idx].onePoint.z = z1;
				NI[idx].onePoint.valide = true;

				NI[idx].twoPoint.x = x2;
				NI[idx].twoPoint.y = y2;
				NI[idx].twoPoint.z = z2;
				NI[idx].twoPoint.valide = true;
			}

			else if (D == 0)
			{
				double SQ = fast_math::sqrt(D);
				double t1 = (-B + SQ) / (2 * A);
				//Первая точка пересечения
				double x1 = L(idx[1]).onePoint.x + t1 * (L(idx[1]).twoPoint.x - L(idx[1]).onePoint.x);
				double y1 = L(idx[1]).onePoint.y + t1 * (L(idx[1]).twoPoint.y - L(idx[1]).onePoint.y);
				double z1 = L(idx[1]).onePoint.z + t1 * (L(idx[1]).twoPoint.z - L(idx[1]).onePoint.z);
					
				NI[idx].onePoint.x = x1;
				NI[idx].onePoint.y = y1;
				NI[idx].onePoint.z = z1;
				NI[idx].onePoint.valide = true;

				NI[idx].twoPoint.valide = false;

			}
			else
			{
				NI[idx].onePoint.valide = false;
				NI[idx].twoPoint.valide = false;
			}
			
		
    });
	NI.synchronize();
	
}


void CPUCalc(std::vector<Line>& Lines, double R, vector<Point>& ScopeCenter, vector<TwoPoints>& NeedIntersections)
{

	for (int j = 0; j < ScopeCenter.size(); j++) 
	{
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

				//TODO:Первая точка пересечения
				double x1 = Lines.at(i).onePoint.x + t1 * (Lines.at(i).twoPoint.x - Lines.at(i).onePoint.x);
				double y1 = Lines.at(i).onePoint.y + t1 * (Lines.at(i).twoPoint.y - Lines.at(i).onePoint.y);
				double z1 = Lines.at(i).onePoint.z + t1 * (Lines.at(i).twoPoint.z - Lines.at(i).onePoint.z);

				//TODO:Вторая точка пересечения
				double x2 = Lines.at(i).onePoint.x + t2 * (Lines.at(i).twoPoint.x - Lines.at(i).onePoint.x);
				double y2 = Lines.at(i).onePoint.y + t2 * (Lines.at(i).twoPoint.y - Lines.at(i).onePoint.y);
				double z2 = Lines.at(i).onePoint.z + t2 * (Lines.at(i).twoPoint.z - Lines.at(i).onePoint.z);


				Point a1;
				a1.x = x1;
				a1.y = y1;
				a1.z = z1;
				
				NeedIntersections[j * Lines.size() + i].onePoint = a1;

				Point a2;
				a2.x = x2;
				a2.y = y2;
				a2.z = z2;
				NeedIntersections[j * Lines.size() + i].onePoint = a2;
			}

			else if (D == 0)
			{
				double t1 = (-B + SQ) / (2 * A);
				//TODO: Первая точка пересечения
				double x1 = Lines.at(i).onePoint.x + t1 * (Lines.at(i).twoPoint.x - Lines.at(i).onePoint.x);
				double y1 = Lines.at(i).onePoint.y + t1 * (Lines.at(i).twoPoint.y - Lines.at(i).onePoint.y);
				double z1 = Lines.at(i).onePoint.z + t1 * (Lines.at(i).twoPoint.z - Lines.at(i).onePoint.z);
				
				Point a1;
				a1.x = x1;
				a1.y = y1;
				a1.z = z1;
				a1.valide = true;

				NeedIntersections[j * Lines.size() + i].onePoint = a1;

				NeedIntersections[j * Lines.size() + i].twoPoint.valide = false;
			
			}
			else
			{
				NeedIntersections[j * Lines.size() + i].onePoint.valide = false;
				NeedIntersections[j * Lines.size() + i].twoPoint.valide = false;
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


	vector <Point> ScopeCenter;
	accelerator defaultDevice(accelerator::default_accelerator);
	accelerator_view defaultView = defaultDevice.default_view;
	wcout << L" Using device : " << defaultDevice.get_description() << endl << endl; //TODO: Название использующейся видеокарты
	
	cout << "Количество шагов по оси X" << endl;
	cin >> StepX;
	
	cout << "Количество шагов по оси Y" << endl;
	cin >> StepY;
	
	const int ArraySize = StepX + StepY + 2; //TODO: Ну тут понятно почему +2, типа количестов линий равно количеству шагов + 1 с каждой стороны
	
	cout << "Введите радиус" << endl;
	cin >> R;
	

	cout << "Колличество сфер " ;
	cin >> NumberOfScope;

	for (int i = 0; i < NumberOfScope; i++)
	{
		//Центры сфер
		Point C;
			
		C.x = 150;
		C.y = 150;
		C.z = 0;

		ScopeCenter.push_back(C);  
	}


	int NN = ScopeCenter.size();

	vector<TwoPoints> NeedIntersections(NN * ArraySize);

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



	NeedIntersections.clear();

	NeedIntersections = vector<TwoPoints>(NN * ArraySize);


	cout << "СPU" << endl;
	unsigned int start_time2 = clock(); //начальное время

	CPUCalc(Lines, R, ScopeCenter, NeedIntersections);

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