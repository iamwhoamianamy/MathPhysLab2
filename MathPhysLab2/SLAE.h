#pragma once
#include <vector>

using namespace std;

// ����� �������� ���� ������� LU - ����������
class SLAE
{
public:
   vector<double> bot_tr;     // ������ ����������� ������� �������
   vector<double> top_tr;     // ������� ����������� ������� �������
   vector<double> di;         // ��������� ������� �������
   vector<int> ind;           // ��������� ������ �����

   int size;                  // ������ ������� �������

   // ��������� ������� �� ������ vec, ��������� � res
   void MatVecMult(const vector<double> &vec, vector<double> &res)
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

   // �������� ������� �������
   void ForwardSolver(vector<double>& res)
   {
      for(int i = 0; i < size; i++)
      {
         int i0 = ind[i + 0], i1 = ind[i + 1];

         double s = 0;

         for(int j = i - (i1 - i0), k = i0; j < i; j++, k++)
            s += res[j] * bot_tr[k];

         res[i] = (res[i] - s) / di[i];
      }
   }

   // �������� ��������� �������
   void BackwardSolver(vector<double>& res)
   {
      for(int i = size - 1; i >= 0; i--)
      {
         int i0 = ind[i + 0], i1 = ind[i + 1];

         double xi = res[i];
         for(int j = i - (i1 - i0), k = i0; j < i; j++, k++)
            res[j] -= xi * top_tr[k];

         res[i] = xi;
      }
   }

   // LU - ���������� ������� �������
   void LUDecomp()
   {
      for(int i = 0; i < size; i++)
      {
         int i0 = ind[i + 0], i1 = ind[i + 1];

         double sd = 0;
         for(int j = i - (i1 - i0), k = i0; j < i; j++, k++)
         {
            double sl = 0, su = 0;
            int j0 = ind[j + 0], j1 = ind[j + 1];
            int kol_i = k - i0, kol_j = j1 - j0;
            int kol_r = kol_i - kol_j, ki = i0, kj = j0;

            if(kol_r > 0)
               ki += kol_r;
            else
               kj -= kol_r;

            for(; ki < k; ki++, kj++)
            {
               sl += bot_tr[ki] * top_tr[kj];
               su += bot_tr[kj] * top_tr[ki];
            }

            bot_tr[k] = bot_tr[k] - sl;
            top_tr[k] = top_tr[k] - su;
            top_tr[k] /= di[j];

            sd += bot_tr[k] * top_tr[k];
         }

         di[i] = di[i] - sd;
      }
   }
}; 