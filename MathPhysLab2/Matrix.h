#pragma once
#include <vector>

using namespace std;

// ����� ������� � ���������� �������
class Matrix
{
public:
   vector<double> bot_tr;     // ������ ����������� ������� �������
   vector<double> top_tr;     // ������� ����������� ������� �������
   vector<double> di;         // ��������� ������� �������
   vector<int> ind;           // ������ ���������� ������ �����
   int size;                  // ������ �������

   // ��������� ������� �� ������ vec, ��������� � res
   void MatVecMult(const vector<double>& vec, vector<double>& res)
   {
      for(int i = 0; i < size; i++)
      {
         res[i] = vec[i] * di[i];

         int i0 = ind[i + 0], i1 = ind[i + 1];

         for(int j = i - (i1 - i0), k = i0; j < i; j++, k++)
         {
            res[i] += vec[j] * bot_tr[k];
            res[j] += vec[i] * top_tr[k];
         }
      }
   }
};
