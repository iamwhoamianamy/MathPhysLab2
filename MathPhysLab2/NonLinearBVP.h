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
   vector<Region> regions;             // Регионы расчетной области

   int n_regions = 0;                  // Количество регионов
   int n_nodes = 0;                    // Общее количество узлов

   SLAE* slae;                         // Система
   Test test;                          // Тестовая информация

   NonLinearBVP()
   {

   }

   // Функция считывания областей из файла FILE_NAME
   // и формирования сетки
   void ReadFormGrid(const string& file_name)
   {
      ifstream fin(file_name);

      fin >> n_regions;
      string s;

      regions.resize(n_regions);

      for(int reg_i = 0; reg_i < n_regions; reg_i++)
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

         r->count = n + 1;
         r->nodes.resize(r->count);

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
            r->first_i = regions[reg_i - 1].count;
      }
      fin.close();
   }
};