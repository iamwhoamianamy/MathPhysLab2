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
      return -1 * (dudx(x, t) * dlambdadx(x) + d2udx2(x, t) * lambda(u_prec(x, t))) + sigma() * dudt(x, t);
   }

   double lambda(const double& u)
   {
      switch(M)
      {
         case 0: return 1;
         case 1: return 3 * u * u + 1;;
      }
   }

   // Производная lambda по x после подстановки x в функцию u
   double dlambdadx(const double& x)
   {
      switch(M)
      {
         case 0: return 0;
         case 1:
         {
            switch(N)
            {
               case(0): return 2.0;
               case(1): return 150 * x;
               case(2): return 12 * x * x * x;
               case(3): return 18 * x * x * x * x * x;
               case(4): return 24 * pow(x, 7);
            }
         }
      }
   }

   double sigma()
   {
      return 1;
   }

   // Точное решение
   double u_prec(const double& x, const double& t)
   {
      switch(N)
      {
         case(0): return 2.0;
         case(1): return t * x + 1;
         case(2): return x * x;
         case(3): return x * x * x;
         case(4): return x * x * x * x;
      };
   }

   // Первая производная точного решения по х
   double dudx(const double& x, const double& t)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return t;
         case(2): return 2 * x;
         case(3): return 3 * x * x;
         case(4): return 4 * x * x * x;
      };
   }

   // Вторая производная точного решения по х
   double d2udx2(const double& x, const double& t)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return 0;
         case(2): return 2;
         case(3): return 6 * x;
         case(4): return 12 * x * x;
      };
   }

   // Производная точного решения по t
   double dudt(const double& x, const double& t)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return x;
         case(2): return 0;
         case(3): return 0;
         case(4): return 0;
      };
   }
};