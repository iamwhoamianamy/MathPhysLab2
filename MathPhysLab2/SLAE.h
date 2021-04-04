#pragma once
#include <vector>
#include "Matrix.h"

using namespace std;

// Класс решателя СЛАУ методом LU - разложения
class SLAE
{
public:
   
   Matrix m;                  // Матрица системы
   vector<double> b;          // Вектор правой части

   // Алгоритм прямого прохода
   void ForwardSolver()
   {
      for(int i = 0; i < m.size; i++)
      {
         int i0 = m.ind[i + 0], i1 = m.ind[i + 1];

         double s = 0;

         for(int j = i - (i1 - i0), k = i0; j < i; j++, k++)
            s += b[j] * m.bot_tr[k];

         b[i] = (b[i] - s) / m.di[i];
      }
   }

   // Алгоритм обратного прохода
   void BackwardSolver()
   {
      for(int i = m.size - 1; i >= 0; i--)
      {
         int i0 = m.ind[i + 0], i1 = m.ind[i + 1];

         double xi = b[i];
         for(int j = i - (i1 - i0), k = i0; j < i; j++, k++)
            b[j] -= xi * m.top_tr[k];

         b[i] = xi;
      }
   }

   // LU - разложение матрицы системы
   void LUDecomp()
   {
      for(int i = 0; i < m.size; i++)
      {
         int i0 = m.ind[i + 0], i1 = m.ind[i + 1];

         double sd = 0;
         for(int j = i - (i1 - i0), k = i0; j < i; j++, k++)
         {
            double sl = 0, su = 0;
            int j0 = m.ind[j + 0], j1 = m.ind[j + 1];
            int kol_i = k - i0, kol_j = j1 - j0;
            int kol_r = kol_i - kol_j, ki = i0, kj = j0;

            if(kol_r > 0)
               ki += kol_r;
            else
               kj -= kol_r;

            for(; ki < k; ki++, kj++)
            {
               sl += m.bot_tr[ki] * m.top_tr[kj];
               su += m.bot_tr[kj] * m.top_tr[ki];
            }

            m.bot_tr[k] = m.bot_tr[k] - sl;
            m.top_tr[k] = m.top_tr[k] - su;
            m.top_tr[k] /= m.di[j];

            sd += m.bot_tr[k] * m.top_tr[k];
         }

         m.di[i] = m.di[i] - sd;
      }
   }
}; 