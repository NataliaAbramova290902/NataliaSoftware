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
	bool valide = true;   //TODO: ���� ���������� �����, ���� � ������ �� ��� ����� ������ ������� ��������� - �������� �����, ���� ��������� ������
};

//TODO: ��������� ��� ��������� ����� �����������, ���� �� �� ����� ���� ������ ����
struct TwoPoints
{
	Point onePoint;
	Point twoPoint;	
};


struct Line
{
	double x = 0;
	double y = 0;
};

void SimpleSqrtGPUCalc(std::vector<Line>& Lines, int StepX, int StepY, double R, vector<Point>& ScopeCenter, vector<TwoPoints>& NeedIntersections)
{
	int N = Lines.size();
	array_view<const Line, 2 > L(StepX, StepY, Lines);
	array_view<const Point, 1 > SC(ScopeCenter.size(), ScopeCenter);
	array_view<TwoPoints, 2 > NI(StepX, StepY, NeedIntersections);
	NI.discard_data();
	

	parallel_for_each(NI.extent, [=](index<2> idx) restrict(amp)
	{

			double A = ((-R - 0.1) - (R + 0.1)) * ((-R - 0.1) - (R + 0.1));
			double B = 2 * (((-R - 0.1) - (R + 0.1)) * ((R + 0.1) - SC(idx[0]).z));
			double C = ((L(idx).x - SC(idx[0]).x) * (L(idx).x - SC(idx[0]).x)) + ((L(idx).y - SC(idx[0]).y) * (L(idx).y - SC(idx[0]).y)) + (((R + 0.1) - SC(idx[0]).z) *((R + 0.1) - SC(idx[0]).z)) - (R * R);

			double D = (B * B) - (4 * A * C);

			if (D > 0)
			{
				double SQ = fast_math::sqrt(D);
				double t1 = (-B + SQ) / (2 * A);
				double t2 = (-B - SQ) / (2 * A);

				//������ ����� �����������
				double x1 = L(idx).x;
				double y1 = L(idx).y;
				double z1 = (R + 0.1) + t1 * ((-R - 0.1) - (R + 0.1));

				//������ ����� �����������
				double x2 = L(idx).x;
				double y2 = L(idx).y;
				double z2 = (R + 0.1) + t2 * ((-R - 0.1) - (R + 0.1));

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
				//������ ����� �����������
				double x1 = L(idx).x;
				double y1 = L(idx).y;
				double z1 = (R + 0.1) + t1 * ((-R - 0.1) - (R + 0.1));
					
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

			double A = ((-R - 0.1) - (R + 0.1)) * ((-R - 0.1) - (R + 0.1));
			double B = 2 * (((-R - 0.1) - (R + 0.1)) * ((R + 0.1) - ScopeCenter.at(j).z));
			double C = ((R + 0.1) - ScopeCenter.at(j).z) * ((R + 0.1) - ScopeCenter.at(j).z) - (R * R);

			double D = (B * B) - (4 * A * C);
			double SQ = fast_math::sqrt(D);
			if (D > 0)
			{
				double t1 = (-B + SQ) / (2 * A);
				double t2 = (-B - SQ) / (2 * A);

				//TODO:������ ����� �����������
				double x1 = Lines.at(i).x;
				double y1 = Lines.at(i).y;
				double z1 = (R + 0.1) + t1 * ((-R - 0.1) - (R + 0.1));

				//TODO:������ ����� �����������
				double x2 = Lines.at(i).x;
				double y2 = Lines.at(i).y;
				double z2 = (R + 0.1) + t2 * ((-R - 0.1) - (R + 0.1));


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
				//TODO: ������ ����� �����������
				double x1 = Lines.at(i).x;
				double y1 = Lines.at(i).y;
				double z1 = (R + 0.1) + t1 * ((-R - 0.1) - (R + 0.1));
				
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
	int NumberOfScope; //TODO: ����������� ����


	//vector <Point> ScopeCenter;
	accelerator defaultDevice(accelerator::default_accelerator);
	accelerator_view defaultView = defaultDevice.default_view;
	wcout << L" Using device : " << defaultDevice.get_description() << endl << endl; //TODO: �������� �������������� ����������

	cout << "���������� ����� �� ��� X" << endl;
	cin >> StepX;

	cout << "���������� ����� �� ��� Y" << endl;
	cin >> StepY;

	const int ArraySize = (StepX + 1) * (StepY + 1);

	cout << "������� ������" << endl;
	cin >> R;


	cout << "����������� ���� (������������ ��������� " << ArraySize << ")" << endl;
	cin >> NumberOfScope;

	vector<Line> Lines(ArraySize); //��������� �����	
	for (int i = 0; i < StepY + 1; i++)
	{
		for (int j = 0; j < StepX + 1; j++)
		{
			Lines[i * (StepY + 1) + j].x = j * (R + 0.1); //���������� X
			Lines[i * (StepY + 1) + j].y = i * (R + 0.1); //���������� Y 
		}
	}
		vector <Point> ScopeCenter(NumberOfScope); //������ ����	
		for (int i = 0; i < NumberOfScope; i++)
		{
			ScopeCenter[i].x = Lines[i].x;
			ScopeCenter[i].y = Lines[i].y;
			ScopeCenter[i].z = 0;
		}

		int NN = ScopeCenter.size();
		vector<TwoPoints> NeedIntersections(ArraySize);

		cout << "GPU" << std::endl << std::endl;
		unsigned int start_time = clock(); //��������� �����

		SimpleSqrtGPUCalc(Lines, StepX, StepX, R, ScopeCenter, NeedIntersections);

		unsigned int end_time = clock(); // �������� �����
		unsigned int search_time = end_time - start_time; // ������� �����
		cout << search_time << " ms." << endl;



		NeedIntersections.clear();

		NeedIntersections = vector<TwoPoints>(ArraySize);


		cout << "�PU" << endl;
		unsigned int start_time2 = clock(); //��������� �����

		CPUCalc(Lines, R, ScopeCenter, NeedIntersections);

		unsigned int end_time2 = clock(); // �������� �����
		unsigned int search_time2 = end_time2 - start_time2; // ������� �����
		cout << search_time2 << " ms." << endl;

		/*for (int i = 0; i < Lines.size(); i++)
		{
			for (int j = 0; j < Lines.at(i).IntersectionPoints.size(); j++)
			{
				cout << "����� " << i + 1 << ", ����� " << Lines[i].IntersectionPoints[j].Scope + 1 << ", X = " << Lines[i].IntersectionPoints[j].x << ", Y = " << Lines.at(i).IntersectionPoints.at(j).y << ", Z = " << Lines.at(i).IntersectionPoints.at(j).z << "\n";
			}
		}*/

		system("pause");
	}
