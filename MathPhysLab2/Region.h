#pragma once
#include <vector>

using namespace std;

struct Region
{
   vector<double> nodes;         // ���������� �����
   int count;                    // ���������� �����

   int first_i;                  // ������ ������� ���� � ���������� ���������

   int left_bord, right_bord;    // ���������� � ������� ��������

};