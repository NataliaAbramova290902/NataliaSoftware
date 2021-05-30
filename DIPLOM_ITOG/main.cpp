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
	
	array_view<const Capsule, 1 > SC(ScopeCenter.size(), ScopeCenter);
	array_view<TwoPoints, 2> NI((StepX + 1) , (StepY + 1), NeedIntersections);
	NI.discard_data();
	
	for (int s = 0; s < ScopeCenter.size(); s++)
	{

		//Перенос системы координат	

		//Вектор который станет осью OX
		Point Vector;
		Vector.x = SC[s].onePoint.x - SC[s].twoPoint.x;
		Vector.y = SC[s].onePoint.y - SC[s].twoPoint.y;
		Vector.z = SC[s].onePoint.z - SC[s].twoPoint.z;

		//Находим углы поворота
		double SinA = Vector.y / (sqrt((Vector.y * Vector.y) + (Vector.x * Vector.x)));
		double CosA = Vector.x / (sqrt((Vector.y * Vector.y) + (Vector.x * Vector.x)));
		double SinB = Vector.z;
		double CosB = sqrt((Vector.y * Vector.y) + (Vector.x * Vector.x));

		//Новые точки капсулы
		Point SC_One;
		SC_One.x = 0;
		SC_One.y = 0;
		SC_One.z = 0;

		Point SC_Two;
		SC_Two.x = abs(sqrt((Vector.x * Vector.x) + (Vector.y * Vector.y) + (Vector.z * Vector.z)));
		SC_Two.y = 0;
		SC_Two.z = 0;

		//Опряделяем прямоугольник
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

		//Для каждого декселя входящего в квадрат
		concurrency::extent<2> e(NI2y - NI1y, NI2x - NI1x);
		parallel_for_each(e, [=](index<2> idx) restrict(amp)
		{
			//исходные координаты
			int j = (NI1x + idx[1]); 
			int i = (NI1y + idx[0]);

			Point NI_New1;
			NI_New1.x = j;
			NI_New1.y = i;
			NI_New1.z = -100000 * R;

			Point NI_New2;
			NI_New2.x = j;
			NI_New2.y = i;
			NI_New2.z = 100000 * R;

			//Координаты после пересноса системы
			//Первая точка
			
			//Перенос в новое начало координат
			NI_New1.x = NI_New1.x + SC[s].onePoint.x;
			NI_New1.y = NI_New1.y + SC[s].onePoint.y;
			NI_New1.z = NI_New1.z + SC[s].onePoint.z;

			//Поворот вокруг оси OZ
			NI_New1.x = NI_New1.x * CosA + NI_New1.y * (-SinA) + NI_New1.z * 0;
			NI_New1.y = NI_New1.x * SinA + NI_New1.y * CosA + NI_New1.z * 0;
			NI_New1.z = NI_New1.x * 0 + NI_New1.y * 0 + NI_New1.z * 1;

			//Поворот вокруг оси OY
			NI_New1.x = NI_New1.x * CosB + NI_New1.y * 0 + NI_New1.z * (-SinB);
			NI_New1.y = NI_New1.x * 0 + NI_New1.y * 1 + NI_New1.z * 0;
			NI_New1.z = NI_New1.x * SinB + NI_New1.y * 0 + NI_New1.z * CosB;

			//Вторая точка
			//Перенос в новое начало координат
			NI_New2.x = NI_New2.x + SC[s].onePoint.x;
			NI_New2.y = NI_New2.y + SC[s].onePoint.y;
			NI_New2.z = NI_New2.z + SC[s].onePoint.z;

			//Поворот вокруг оси OZ
			NI_New2.x = NI_New2.x * CosA + NI_New2.y * (-SinA) + NI_New2.z * 0;
			NI_New2.y = NI_New2.x * SinA + NI_New2.y * CosA + NI_New2.z * 0;
			NI_New2.z = NI_New2.x * 0 + NI_New2.y * 0 + NI_New2.z * 1;

			//Поворот вокруг оси OY
			NI_New2.x = NI_New2.x * CosB + NI_New2.y * 0 + NI_New2.z * (-SinB);
			NI_New2.y = NI_New2.x * 0 + NI_New2.y * 1 + NI_New2.z * 0;
			NI_New2.z = NI_New2.x * SinB + NI_New2.y * 0 + NI_New2.z * CosB;
				
			double A;
			double B;
			double C;

			if (NI_New1.x < SC_One.x)//первая сфера
			{
				A = (NI_New2.z - NI_New1.z) * (NI_New2.z - NI_New1.z);
				B = 2 * ((NI_New2.z - NI_New1.z) * (NI_New1.z - SC_One.z));
				C = ((NI_New1.x - SC_One.x) * (NI_New1.x - SC_One.x)) + ((NI_New1.y - SC_One.y) * (NI_New1.y - SC_One.y)) + ((NI_New1.z - SC_One.z) * (NI_New1.z - SC_One.z)) - (R * R);
			}
			else if (NI_New1.x > SC_Two.x)
			{
				A = (NI_New2.z - NI_New1.z) * (NI_New2.z - NI_New1.z);
				B = 2 * ((NI_New2.z - NI_New1.z) * (NI_New1.z - SC_Two.z));
				C = ((NI_New1.x - SC_Two.x) * (NI_New1.x - SC_Two.x)) + ((NI_New1.y - SC_Two.y) * (NI_New1.y - SC_Two.y)) + ((NI_New1.z - SC_Two.z) * (NI_New1.z - SC_Two.z)) - (R * R);
			}

			else
			{
				A = (NI_New2.z - NI_New1.z) * (NI_New2.z - NI_New1.z);
				B = 2 * ((NI_New2.z - NI_New1.z) * (NI_New1.z - SC_One.z));
				C = ((NI_New1.y - SC_One.y) * (NI_New1.y - SC_One.y)) + ((NI_New1.z - SC_One.z) * (NI_New1.z - SC_One.z)) - (R * R);
			}

			double D = (B * B) - 4 * A * C;

			if (D > 0)
			{
				double SQ = fast_math::sqrt(D);
				double t1 = (-B + SQ) / (2 * A);
				double t2 = (-B - SQ) / (2 * A);

				//Первая точка пересечения
				NI_New1.z = NI_New1.z + t1 * (NI_New2.z - NI_New1.z);

				//Вторая точка пересечения
				NI_New2.z = NI_New1.z + t2 * (NI_New2.z - NI_New1.z);

				//Преобразовать координаты обратно

				//Первая точка

				//Поворот вокруг оси OY
				NI_New1.z = NI_New1.x * (-SinB) + NI_New1.y * 0 + NI_New1.z * CosB;

				//Поворот вокруг оси OZ
				NI_New1.z = NI_New1.x * 0 + NI_New1.y * 0 + NI_New1.z * 1;

				//Перенос в новое начало координат
				NI_New1.z = NI_New1.z - SC[s].onePoint.z;

				//Вторая точка
				//Поворот вокруг оси OY
				NI_New2.z = NI_New2.x * (-SinB) + NI_New2.y * 0 + NI_New2.z * CosB;

				//Поворот вокруг оси OZ
				NI_New2.z = NI_New2.x * 0 + NI_New2.y * 0 + NI_New2.z * 1;

				//Перенос в новое начало координат
				NI_New2.z = NI_New2.z - SC[s].onePoint.z;

				if (NI_New1.z < NI_New2.z)
				{
					if(NI[i][j].onePointZ > NI_New1.z)
						NI[i][j].onePointZ = NI_New1.z;

					if(NI[i][j].twoPointZ < NI_New2.z)
						NI[i][j].twoPointZ = NI_New2.z;
				}
				else
				{
					if (NI[i][j].onePointZ > NI_New2.z)
						NI[i][j].onePointZ = NI_New2.z;

					if (NI[i][j].twoPointZ < NI_New1.z)
						NI[i][j].twoPointZ = NI_New1.z;
				}
			}

			else if (D == 0)
			{
				double SQ = fast_math::sqrt(D);
				double t1 = (-B + SQ) / (2 * A);
				//Первая точка пересечения
				NI_New1.z = NI_New1.z + t1 * (NI_New2.z - NI_New1.z);

				NI_New2.z = NI_New1.z + t1 * (NI_New2.z - NI_New1.z);

				//Преобразовать координаты обратно

				//Первая точка

				//Поворот вокруг оси OY
				NI_New1.z = NI_New1.x * (-SinB) + NI_New1.y * 0 + NI_New1.z * CosB;

				//Поворот вокруг оси OZ
				NI_New1.z = NI_New1.x * 0 + NI_New1.y * 0 + NI_New1.z * 1;

				//Перенос в новое начало координат
				NI_New1.z = NI_New1.z - SC[s].onePoint.z;

				if (NI[i][j].onePointZ > NI_New1.z)
					NI[i][j].onePointZ = NI_New1.z;

				if (NI[i][j].twoPointZ < NI_New1.z)
					NI[i][j].twoPointZ = NI_New1.z;

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

	cout << "Количество шагов" << endl;
	cin >> StepX;

	StepY = StepX;

	const int ArraySize = (StepX + 1) * (StepY + 1);

	cout << "Введите радиус инструмента" << endl;
	cin >> R;

	cout << "Колличество капсул" << endl;
	cin >> NumberOfScope;

	vector <Capsule> ScopeCenter(NumberOfScope); //Точки капсул

	for (int i = 0; i < NumberOfScope; i++)
	{
		ScopeCenter[i] = generateCapsule(R, 5*R, StepX);

	}

	vector<TwoPoints> NeedIntersections(ArraySize);

	cout << "GPU_v3" << std::endl;
	unsigned int start_time = clock(); //начальное время

	GPUCalc_3(StepX, StepX, R, ScopeCenter, NeedIntersections);

	unsigned int end_time = clock(); // конечное время
	unsigned int search_time = end_time - start_time; // искомое время
	cout << search_time << " ms." << endl;

	NeedIntersections.clear();

	system("pause");
}
