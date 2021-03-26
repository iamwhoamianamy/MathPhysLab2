#pragma once
#include <vector>

using namespace std;

struct Region
{
   vector<double> nodes;         // Координаты узлов
   int nodes_count;              // Количество узлов
   int elems_count;              // Количество конечных элементов

   int first_i;                  // Индекс первого узла в глобальной нумерации

   int left_bord, right_bord;    // Информация о краевых условиях

};