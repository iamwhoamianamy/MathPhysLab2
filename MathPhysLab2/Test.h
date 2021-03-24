#pragma once
#pragma once
using namespace std;

class Test
{
public:

   int N;

   Test(const int& t_N) : N(t_N) {};

   Test() : N(0) {};

   double f(const double& x, const double& y)
   {
      switch(N)
      {
         case(0): return (0) * lambda() + u(x, y) * gamma();
         case(1): return (0) * lambda() + u(x, y) * gamma();
         case(2): return (-4) * lambda() + u(x, y) * gamma();
         case(3): return (-6 * x - 6 * y) * lambda() + u(x, y) * gamma();
         case(4): return (-12 * x * x - 12 * y * y) * lambda() + u(x, y) * gamma();
      };
   }

   double lambda()
   {
      return 1;
   }

   double gamma()
   {
      return 1;
   }

   double u(const double& x, const double& y)
   {
      switch(N)
      {
         case(0): return 2.0;
         case(1): return x + y;
         case(2): return x * x + y * y;
         case(3): return x * x * x + y * y * y;
         case(4): return x * x * x * x + y * y * y * y;
      };
   }
};