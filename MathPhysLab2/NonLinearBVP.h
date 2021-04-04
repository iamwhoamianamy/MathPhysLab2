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
   vector<double> q_1;        // Приближение функции u на предыдущей итерации
   vector<double> q;          // Приближение функции u на текущей итерации

   Matrix G, M;               // Матрицы жесткости и массы
   SLAE slae;                 // СЛАУ
   Test test;                 // Тестовая информация

   vector<double> time_grid;  // Сетка по времени

   vector<double> vec_1;      // Вспомогательный вектор для метода
                              // простой итерации

   vector<double> vec_2;      // Вспомогательный вектор для явной схемы
                              // по времени

   // Вспомогательная матрица для построения
   // матрицы массы конечного элемента
   vector<vector<int>> C;

   NonLinearBVP()
   {
      C = {
         {4, 2, -1},
         {2, 16, 2},
         {-1, 2, 4}
      };
   }

   // Функция считывания областей из файла file_name
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

      vec_1.resize(nodes_count);
   }

   // Функция считывания и формирования сетки по времени
   // из файла file_name
   void ReadFormtimeGrid(const string& file_name)
   {
      ifstream fin(file_name);

      double t0, tn;

      fin >> t0;
      fin >> tn;

      // Генерация координат узлов по X
      int n;
      double h, q;

      fin >> q >> n;

      h = tn - t0;

      if(q != 1)
         h *= (1 - q) / (1 - pow(q, n));
      else
         h /= n;

      time_grid.resize(n + 1);
      time_grid[0] = t0;

      for(int i = 0; i < n; i++)
         time_grid[i + 1] = time_grid[i] + h * pow(q, i);

      fin.close();
   }

   // Функция инициализации матрицы
   void InitMatrix(Matrix& m)
   {
      m.size = 0;
      m.bot_tr = vector<double>(elems_count * 3);
      m.top_tr = vector<double>(elems_count * 3);

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         m.size += regions[reg_i].elems_count * 2 + 1;

         if(reg_i + 1 < regions_count && grid[reg_i * 2 + 1] == grid[(reg_i + 1) * 2])
            m.size--;
      }

      m.di = vector<double>(m.size);
      m.ind = vector<int>(m.size + 1);
   }

   // Функция формирования портрета глобальной матрицы
   void FormPortrait()
   {
      slae.m.ind[0] = 0;
      slae.m.ind[1] = 0;
      slae.m.ind[2] = 1;
      slae.m.ind[slae.m.size] = slae.m.top_tr.size();

      int global_i = 3;
      int val = 3;

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         Region* r = &regions[reg_i];

         for(int elem_i = reg_i == 0 ? 1 : 0; elem_i < r->elems_count; elem_i++, global_i += 2, val += 2)
         {
            slae.m.ind[global_i] = val++;
            slae.m.ind[global_i + 1] = val;
         }

         if(reg_i + 1 < regions_count && grid[reg_i * 2 + 1] != grid[(reg_i + 1) * 2])
            slae.m.ind[global_i++] = val;
      }
   }

   // Функция заполнения матриц жесткости и массы
   void FillMatrices()
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

            // Заполнение диагонали матрицы жесткости
            G.di[to_add_i_di++] += (lambda[0] * 37 / 30 + lambda[1] * 1.2 - lambda[2] * 0.1) / h;
            G.di[to_add_i_di++] += (lambda[0] * 1.6 + lambda[1] * 32 / 15 + lambda[2] * 1.6) / h;
            G.di[to_add_i_di]   += (-lambda[0] * 0.1 + lambda[1] * 1.2 + lambda[2] * 37 / 30) / h;

            // Заполнение нижнего треугольника матрицы жесткости
            G.bot_tr[to_add_i_tr++] += (-lambda[0] * 22 / 15 - lambda[1] * 16 / 15 - lambda[2] * 2 / 15) / h;
            G.bot_tr[to_add_i_tr++] += (lambda[0] * 7 / 30 - lambda[1] * 2 / 15 + lambda[2] * 7 / 30) / h;
            G.bot_tr[to_add_i_tr++] += (-lambda[0] * 2 / 15 - lambda[1] * 16 / 15 - lambda[2] * 22 / 15) / h;

            to_add_i_di -= 2;
            to_add_i_tr -= 3;

            // Заполнение диагонали матрицы массы
            M.di[to_add_i_di++] += test.gamma() * h / 30.0 * C[0][0];
            M.di[to_add_i_di++] += test.gamma() * h / 30.0 * C[1][1];
            M.di[to_add_i_di]   += test.gamma() * h / 30.0 * C[2][2];

            // Заполнение нижнего треугольника массы
            M.bot_tr[to_add_i_tr++] += test.gamma() * h / 30.0 * C[1][0];
            M.bot_tr[to_add_i_tr++] += test.gamma() * h / 30.0 * C[2][0];
            M.bot_tr[to_add_i_tr++] += test.gamma() * h / 30.0 * C[2][1];

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
      G.top_tr = G.bot_tr;
      M.top_tr = M.bot_tr;
   }
   
   // Функция сбора глобальной матрицы системы
   void GlobalMatrixAssemble()
   {
      int tr_size = slae.m.bot_tr.size();
      for(int i = 0; i < tr_size; i++)
         slae.m.bot_tr[i] = G.bot_tr[i] + M.bot_tr[i];

      slae.m.top_tr = slae.m.bot_tr;

      for(int i = 0; i < slae.m.size; i++)
         slae.m.di[i] = G.di[i] + M.di[i];
   }

   // Функция учета краевых условий
   void AccountBound()
   {
      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         Region* r = &regions[reg_i];

         if(r->left_bord == 1)
         {
            slae.m.di[r->first_i] = big_num;
            slae.b[r->first_i] = big_num * test.u_prec(r->nodes[0]);
         }

         if(r->right_bord == 1)
         {
            slae.m.di[r->first_i + r->nodes_count - 1] = big_num;
            slae.b[r->first_i + r->nodes_count - 1] = big_num * test.u_prec(r->nodes[r->nodes_count - 1]);
         }
      }
   }

   // Функция расчета нормы вектора
   static double CalcNorm(const vector<double>& vec)
   {
      double res = 0;
      int size = vec.size();

      for(int i = 0; i < size; i++)
         res += vec[i] * vec[i];
      
      return sqrt(res);
   }

   // Линеаризация методом простой итераций, вывод итераций в файл file_name
   void SimpleIterations(const double& eps, const double& delta, const int& max_iter, ofstream& fout)
   {
      int n_iter = 0;
      double eps_residual = 1;
      double delta_residual = 1;

      do
      {
         // Инициализация матриц
         InitMatrix(slae.m);
         InitMatrix(G);
         InitMatrix(M);
         slae.b = vector<double>(slae.m.size);

         // Сборка матриц жесткости и массы
         FormPortrait();
         FillMatrices();

         // Сборка матрицы системы
         GlobalMatrixAssemble();
         AccountBound();

         /*vec_1 = vector<double>(nodes_count);
         slae.m.MatVecMult(q, vec_1);

         for(int i = 0; i < nodes_count; i++)
            vec_1[i] -= slae.b[i];

         eps_residual = CalcNorm(vec_1) / CalcNorm(slae.b) / big_num;

         if(eps_residual < eps)
         {
            fout << endl << "Non-linear iteration " << n_iter << endl;
            PrintSolution(fout);
            break;
         }
         else*/
         {
            slae.LUDecomp();
            slae.ForwardSolver();
            slae.BackwardSolver();

            for(int i = 0; i < nodes_count; i++)
               vec_1[i] = slae.b[i] - q[i];

            delta_residual = CalcNorm(vec_1) / CalcNorm(slae.b);

            if(delta_residual < delta)
            {
               break;
            }
            q = slae.b;

            fout << endl << "Non-linear iteration " << n_iter << endl;
            PrintSolution(fout);
            n_iter++;
         }
      } while(0 < n_iter && n_iter <  max_iter);
   }

   // Явная разностная схема по времени, вывод итераций в поток fout
   void ExplicitScheme(const string& file_name)
   {
      ofstream fout(file_name);

      // Расчет вектора приближения при t=0
      q_1.resize(nodes_count);
      vec_2.resize(nodes_count);

      int to_add_i = 0;
      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         Region* r = &regions[reg_i];

         for(int elem_i = 0; elem_i < r->elems_count; elem_i++)
         {
            q_1[to_add_i++] = test.u_prec(r->nodes[elem_i * 2]);
            q_1[to_add_i++] = test.u_prec(r->nodes[elem_i * 2 + 1]);
            q_1[to_add_i] = test.u_prec(r->nodes[elem_i * 2 + 2]);
         }
      }

      for(int time_i = 1; time_i < time_grid.size(); time_i++)
      {
         // Решение нелинейности
         SimpleIterations(1e-14, 1e-14, 100, fout);

         // Инициализация матриц
         InitMatrix(slae.m);
         InitMatrix(G);
         InitMatrix(M);
         slae.b = vector<double>(slae.m.size);

         // Сборка матриц жесткости и массы
         FormPortrait();
         FillMatrices();

         // Сборка матрицы системы

         // Временной щаг
         double dt = time_grid[time_i] - time_grid[time_i - 1];

         // Рачет матриц системы
         int tr_size = slae.m.bot_tr.size();
         for(int i = 0; i < tr_size; i++)
            slae.m.bot_tr[i] = G.bot_tr[i] + M.bot_tr[i] / dt;

         slae.m.top_tr = slae.m.bot_tr;
         for(int i = 0; i < slae.m.size; i++)
            slae.m.di[i] = G.di[i] + M.di[i] / dt;

         // Расчет вектора правой части
         M.MatVecMult(q_1, vec_2);

         for(int i = 0; i < slae.m.size; i++)
            slae.b[i] += vec_2[i] / dt;

         AccountBound();

         slae.LUDecomp();
         slae.ForwardSolver();
         slae.BackwardSolver();
         q_1 = q;
         q = slae.b;

         fout << endl << "Time iteration " << time_i<< endl;
         PrintSolution(fout);
      }


      fout.close();
   }

   // Вывод решения на основе вектора текущего приближения в поток fout
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
            double help_1 = q[elem_beg_i];
            fout << setw(15) << help_1;
            double help_2 = test.u_prec(r->nodes[node_i]);
            fout << setw(15) << help_2;
            fout << setw(14) << scientific << abs(help_1 - help_2) << fixed;

            fout << setw(4) << elem_beg_i << " ";


            if(node_i == 0 || node_i == r->nodes_count - 1)
               fout << "  border";
            else
               fout << "  inner";

            fout << endl;

            norm_u += help_2 * help_2;
            norm += abs(help_1 - help_2) * abs(help_1 - help_2);
         }
      }
      fout << "||u-u*||/||u*|| = " << scientific << sqrt(norm) / sqrt(norm_u) << endl;
      fout << "||u-u*|| = " << scientific << sqrt(norm) << endl;
   }
};