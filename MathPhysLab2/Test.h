#pragma once
#pragma once
using namespace std;

class Test
{
public:

   int N, M;

   Test(const int& t_N, const int& t_M) : N(t_N), M(t_M) {};

   Test() : N(0), M(0) {};

   double f(const double& x, const double& t)
   {
      return -1 * (dudx(x, t) * dlambdadx(x, t) + d2udx2(x, t) * lambda(u_prec(x, t), x, t)) + sigma() * dudt(x, t);
   }

   double sigma()
   {
      return 1;
   }

   double lambda(const double& u, const double& x, const double& t)
   {
      switch(M)
      {
         case 0: return 1;
         //case 1: return u * u + 1;
         case 1: return x * x + 1;
         //case 1: return u + 1;
         //case 1: return u + t * x;
      }
   }

   // Производная lambda по x после подстановки x в функцию u
   double dlambdadx(const double& x, const double& t)
   {
      switch(M)
      {
         case 0: return 0;
         case 1:
         {
            switch(N)
            {
               case(0): return 0;
               //case(1): return 1;
               case(1): return 2 * x;
               //case(1): return 1 + t;
               //case(2): return t;
               case(2): return 2 * x * t * t;
               //case(3): return 2 * x;
               case(3): return 4 * x * x * x;
               case(4): return 3 * x * x;
               //case(4): return 2 * x + 1;
               case(5): return 4 * x * x * x;
               case(6): return 0;
               //case(7): return 2 * x * t * t;
               case(7): return 0;
               case(8): return  2 * x * t;
            }
         }
      }
   }

   // Точное решение
   double u_prec(const double& x, const double& t)
   {
      switch(N)
      {
         case(0): return 2.0;
         case(1): return x;
         case(2): return x * t;
         case(3): return x * x;
         case(4): return x * x * x;
         case(5): return x * x * x * x;
         //case(5): return x * x * x * x * t;
         case(6): return t;
         case(7): return t * t;
         case(8): return x * x * t;
      };
   }

   // Первая производная точного решения по х
   double dudx(const double& x, const double& t)
   {
      switch(N)
      {
         case(1): return 1;
         case(2): return t;
         case(3): return 2 * x;
         case(4): return 3 * x * x;
         case(5): return 4 * x * x * x;
         case(8): return 2 * x * t;
         default: return 0;
      };
   }

   // Вторая производная точного решения по х
   double d2udx2(const double& x, const double& t)
   {
      switch(N)
      {
         case(3): return 2;
         case(4): return 6 * x;
         case(5): return 12 * x * x;
         case(8): return  2 * t;
         default: return 0;
      };
   }

   // Производная точного решения по t
   double dudt(const double& x, const double& t)
   {
      switch(N)
      {
         case(2): return x;
         case(6): return 1;
         //case(5): return x * x * x * x;
         case(7): return 2*t;
         case(8): return x * x;
         default: return 0;
      };
   }
};