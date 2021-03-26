#pragma once
#include <vector>
#include <fstream>
#include "Region.h"
#include "SLAE.h"
#include "Test.h"

using namespace std;

class NonLinearBVP
{
public:
   vector<Region> regions;    // Регионы расчетной области

   int regions_count = 0;         // Количество регионов
   int nodes_count = 0;           // Общее количество узлов
   int elems_count = 0;           // Общее количество конечных элементов

   SLAE slae;                // СЛАУ
   Test test;                 // Тестовая информация

   // Вспомогательные матрицы для построения матриц 
   // жесткости и массы конечного элемента
   vector<vector<int>> G, C;

   NonLinearBVP()
   {
      G = {
         {7, -8, 1},
         {-8, 16, -8},
         {1, -8, 7}
      };

      C = {
         {4, 2, -1},
         {2, 16, 2},
         {-1, 2, 4}
      };

   }

   // Функция считывания областей из файла FILE_NAME
   // и формирования сетки
   void ReadFormGrid(const string& file_name)
   {
      ifstream fin(file_name);

      fin >> regions_count;
      string s;

      regions.resize(regions_count);

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         fin >> s;
         Region* r = &regions[reg_i];

         double left, right;

         // Считывание границы области
         fin >> left;
         fin >> right;

         // Генерация координат узлов по X
         int n;
         double h, q;

         fin >> q >> n;

         r->nodes_count = n + 1;
         r->nodes.resize(r->nodes_count);

         h = right - left;

         if(q != 1)
            h *= (1 - q) / (1 - pow(q, n));
         else
            h /= n;

         r->nodes[0] = left;

         for(int i = 0; i < n; i++)
            r->nodes[i + 1] = r->nodes[i] + h * pow(q, i);

         // Считывание информации о краевых условиях
         fin >> r->left_bord;
         fin >> r->right_bord;

         if(reg_i > 0)
            r->first_i = regions[reg_i - 1].nodes_count;

         nodes_count += r->nodes_count;

         r->elems_count = r->nodes_count / 2;

         elems_count += r->elems_count;
      }
      fin.close();
   }

   // Функция формирования портрета глобальной матрицы
   void FormPortrait()
   {
      slae.bot_tr.resize(elems_count * 3);
      slae.top_tr.resize(elems_count * 3);
      slae.size = elems_count * 2 + 1;
      slae.di.resize(slae.size);
      slae.ind.resize(slae.size + 1);

      slae.ind[0] = 0;
      slae.ind[1] = 0;
      slae.ind[2] = 1;
      slae.ind[slae.size] = slae.top_tr.size();

      int global_i = 3;

      for(int elem_i = 1; elem_i < elems_count; elem_i++, global_i += 2)
      {
         slae.ind[global_i] = (elem_i) * 3;
         slae.ind[global_i + 1] = slae.ind[global_i] + 1;
      }
   }

   // Функция заполнения матрицы системы
   void FillMatrix()
   {

   }
};