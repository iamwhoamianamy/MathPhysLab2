#pragma once
#include <vector>
#include <fstream>
#include <iomanip>
#include "Region.h"
#include "SLAE.h"
#include "Test.h"

using namespace std;

class NonLinearBVP
{
public:
   vector<Region> regions;    // Регионы расчетной области
   vector<double> grid;       // Исходная сетка

   int regions_count = 0;     // Количество регионов
   int nodes_count = 0;       // Общее количество узлов
   int elems_count = 0;       // Общее количество конечных элементов

   double big_num = 1E+20;    // Большое число для учета первого краевого условия
   vector<double> q;         // Приближение функции u

   SLAE slae;                 // СЛАУ
   Test test;                 // Тестовая информация

   vector<double> t;         // Вспомогательный вектор для метода
                              // простой итерации

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
      regions.resize(regions_count);
      grid.resize(regions_count * 2);

      for(int grid_i = 0; grid_i < regions_count * 2; grid_i++)
         fin >> grid[grid_i];

      string s;
      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         fin >> s;
         Region* r = &regions[reg_i];

         int left, right;

         // Считывание границы области
         fin >> left;
         fin >> right;

         // Генерация координат узлов по X
         int n;
         double h, q;

         fin >> q >> n;

         r->nodes_count = n + 1;
         r->nodes.resize(r->nodes_count);

         h = grid[right] - grid[left];

         if(q != 1)
            h *= (1 - q) / (1 - pow(q, n));
         else
            h /= n;

         r->nodes[0] = grid[left];

         for(int i = 0; i < n; i++)
            r->nodes[i + 1] = r->nodes[i] + h * pow(q, i);

         // Считывание информации о краевых условиях
         fin >> r->left_bord;
         fin >> r->right_bord;

         if(reg_i > 0)
            if(grid[reg_i * 2] == grid[(reg_i - 1) * 2 + 1])
               r->first_i = regions[reg_i - 1].nodes_count - 1;
            else
               r->first_i = regions[reg_i - 1].nodes_count;

         nodes_count += r->nodes_count;

         r->elems_count = r->nodes_count / 2;

         elems_count += r->elems_count;
      }
      fin.close();

      q.resize(nodes_count, 0);

      /*for(int i = 0; i < nodes_count; i++)
         q[i] = pow(11 + 0.25 * (i), 2) ;*/

      t.resize(nodes_count);
   }

   // Функция инициализации матриц системы
   void InitMatrix()
   {
      slae.bot_tr = vector<double>(elems_count * 3);
      slae.top_tr = vector<double>(elems_count * 3);
      slae.size = 0;

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         slae.size += regions[reg_i].elems_count * 2 + 1;

         if(reg_i + 1 < regions_count && grid[reg_i * 2 + 1] == grid[(reg_i + 1) * 2])
            slae.size--;
      }

      slae.di = vector<double>(slae.size);
      slae.b = vector<double>(slae.size);
      slae.ind = vector<int>(slae.size + 1);
   }

   // Функция формирования портрета глобальной матрицы
   void FormPortrait()
   {
      slae.ind[0] = 0;
      slae.ind[1] = 0;
      slae.ind[2] = 1;
      slae.ind[slae.size] = slae.top_tr.size();

      int global_i = 3;
      int val = 3;

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         Region* r = &regions[reg_i];

         for(int elem_i = reg_i == 0 ? 1 : 0; elem_i < r->elems_count; elem_i++, global_i += 2, val += 2)
         {
            slae.ind[global_i] = val++;
            slae.ind[global_i + 1] = val;
         }

         if(reg_i + 1 < regions_count && grid[reg_i * 2 + 1] != grid[(reg_i + 1) * 2])
            slae.ind[global_i++] = val;
      }
   }

   // Функция заполнения матрицы системы
   void FillMatrix()
   {
      // Индекс очередного элемента в треугольнике матрицы
      int to_add_i_tr = 0;
      // Индекс очередного элемента на диагонали матрицы
      int to_add_i_di = 0;

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         Region* r = &regions[reg_i];

         for(int elem_i = 0; elem_i < r->elems_count; elem_i++)
         {
            // Индекс первого узла элемента
            int elem_beg_i = elem_i * 2;

            // Координаты узлов
            vector<double> x_elem = { r->nodes[elem_beg_i], r->nodes[elem_beg_i + 1], r->nodes[elem_beg_i + 2] };
            double h = x_elem[2] - x_elem[0];

            // Индекс первого узла элемента в глобальной нумерации
            int elem_beg_i_glob = elem_i * 2 + r->first_i;

            vector<double> q_elem = { q[elem_beg_i_glob], q[elem_beg_i_glob + 1], q[elem_beg_i_glob + 2] };

            vector<double> lambda = { test.lambda(q_elem[0]),
                                      test.lambda(q_elem[1]),
                                      test.lambda(q_elem[2]) };

            // Заполнение диагонали матрицы
            slae.di[to_add_i_di++] += (lambda[0] * 37 / 30 + lambda[1] * 1.2 - lambda[2] * 0.1) / h + 
                                       test.gamma() * h / 30.0 * C[0][0];
            slae.di[to_add_i_di++] += (lambda[0] * 1.6 + lambda[1] * 32 / 15 + lambda[2] * 1.6) / h +
                                       test.gamma() * h / 30.0 * C[1][1];
            slae.di[to_add_i_di] +=   (-lambda[0] * 0.1 + lambda[1] * 1.2 + lambda[2] * 37 / 30) / h +
                                       test.gamma() * h / 30.0 * C[2][2];

            // Заполнение нижнего треугольника матрицы
            slae.bot_tr[to_add_i_tr++] += (-lambda[0] * 22 / 15 - lambda[1] * 16 / 15 - lambda[2] * 2 / 15) / h +
                                           test.gamma() * h / 30.0 * C[1][0];
            slae.bot_tr[to_add_i_tr++] += (lambda[0] * 7 / 30 - lambda[1] * 2 / 15 + lambda[2] * 7 / 30) / h +
                                           test.gamma() * h / 30.0 * C[2][0];
            slae.bot_tr[to_add_i_tr++] += (-lambda[0] * 2 / 15 - lambda[1] * 16 / 15 - lambda[2] * 22 / 15) / h +
                                           test.gamma() * h / 30.0 * C[2][1];

            vector<double> f_elem = { test.f(x_elem[0]),
                                      test.f(x_elem[1]),
                                      test.f(x_elem[2]) };

            // Заполнение вектора правой части
            slae.b[to_add_i_di - 2] += h / 30.0 * (C[0][0] * f_elem[0] +
                                                   C[0][1] * f_elem[1] +
                                                   C[0][2] * f_elem[2]);
            slae.b[to_add_i_di - 1] += h / 30.0 * (C[1][0] * f_elem[0] +
                                                   C[1][1] * f_elem[1] +
                                                   C[1][2] * f_elem[2]);
            slae.b[to_add_i_di - 0] += h / 30.0 * (C[2][0] * f_elem[0] +
                                                   C[2][1] * f_elem[1] +
                                                   C[2][2] * f_elem[2]);
         }
      }
      slae.top_tr = slae.bot_tr;
   }

   // Функция учета краевых условий
   void AccountBound()
   {
      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         Region* r = &regions[reg_i];

         if(r->left_bord == 1)
         {
            slae.di[r->first_i] = big_num;
            slae.b[r->first_i] = big_num * test.u_prec(r->nodes[0]);
         }

         if(r->right_bord == 1)
         {
            slae.di[r->first_i + r->nodes_count - 1] = big_num;
            slae.b[r->first_i + r->nodes_count - 1] = big_num * test.u_prec(r->nodes[r->nodes_count - 1]);
         }
      }
   }

   double CalcNorm(const vector<double>& vec)
   {
      double res = 0;
      int size = vec.size();

      for(int i = 0; i < size; i++)
         res += vec[i] * vec[i];
      
      return sqrt(res);
   }

   // Линеаризация методом простой итераций, вывод итераций в файл file_name
   void SimpleIterations(const double& eps, const double& delta, const int& max_iter, const string& file_name)
   {
      ofstream fout(file_name);

      int n_iter = 0;
      double eps_residual = 1;
      double delta_residual = 1;

      do
      {
         InitMatrix();
         FormPortrait();
         FillMatrix();
         AccountBound();

         if(n_iter > 0)
         {
            slae.MatVecMult(q, t);

            for(int i = 0; i < nodes_count; i++)
               t[i] -= slae.b[i];

            eps_residual = CalcNorm(t) / CalcNorm(slae.b) * big_num;
         }

         slae.LUDecomp();
         slae.ForwardSolver();
         slae.BackwardSolver();

         for(int i = 0; i < nodes_count; i++)
            t[i] = slae.b[i] - q[i];

         delta_residual = CalcNorm(t) / CalcNorm(slae.b);

         fout << endl << "Iteration = " << n_iter << endl;
         PrintSolution(fout);

         q = slae.b;
         n_iter++;
      } while(0 < n_iter && n_iter <  max_iter && eps_residual > eps && delta_residual > delta);

      fout.close();
   }

   // Вывод решения в файл FILE_NAME
   void PrintSolution(const string& file_name)
   {
      ofstream fout(file_name);
      double norm = 0., norm_u = 0.;

      fout << "   x              calc           prec      dif            N  location" << endl << fixed;

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         Region* r = &regions[reg_i];

         for(int node_i = 0; node_i < r->nodes_count; node_i++)
         {
            // Индекс узла в глобальной нумерации
            int elem_beg_i = node_i + r->first_i;

            fout << setw(11) << r->nodes[node_i];
            double t = slae.b[elem_beg_i];
            fout << setw(15) << t;
            double tt = test.u_prec(r->nodes[node_i]);
            fout << setw(15) << tt;
            fout << setw(14) << scientific << abs(t - tt) << fixed;

            fout << setw(4) << elem_beg_i << " ";

           
           if(node_i == 0 || node_i == r->nodes_count - 1)
               fout << "  border";
            else
               fout << "  inner";

            fout << endl;

            norm_u += tt * tt;
            norm += abs(t - tt) * abs(t - tt);
         }
      }
      fout << "||u-u*||/||u*|| = " << scientific << sqrt(norm) / sqrt(norm_u) << endl;
      fout << "||u-u*|| = " << scientific << sqrt(norm) << endl;
      fout.close();
   }

   // Вывод решения в файл поток fout
   void PrintSolution(ofstream& fout)
   {
      double norm = 0., norm_u = 0.;

      fout << "   x              calc           prec      dif            N  location" << endl << fixed;

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         Region* r = &regions[reg_i];

         for(int node_i = 0; node_i < r->nodes_count; node_i++)
         {
            // Индекс узла в глобальной нумерации
            int elem_beg_i = node_i + r->first_i;

            fout << setw(11) << r->nodes[node_i];
            double t = slae.b[elem_beg_i];
            fout << setw(15) << t;
            double tt = test.u_prec(r->nodes[node_i]);
            fout << setw(15) << tt;
            fout << setw(14) << scientific << abs(t - tt) << fixed;

            fout << setw(4) << elem_beg_i << " ";


            if(node_i == 0 || node_i == r->nodes_count - 1)
               fout << "  border";
            else
               fout << "  inner";

            fout << endl;

            norm_u += tt * tt;
            norm += abs(t - tt) * abs(t - tt);
         }
      }
      fout << "||u-u*||/||u*|| = " << scientific << sqrt(norm) / sqrt(norm_u) << endl;
      fout << "||u-u*|| = " << scientific << sqrt(norm) << endl;
   }
};