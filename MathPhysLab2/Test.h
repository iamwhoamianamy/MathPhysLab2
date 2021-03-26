#pragma once
#pragma once
using namespace std;

class Test
{
public:

   int N;

   Test(const int& t_N) : N(t_N) {};

   Test() : N(0) {};

   double f(const double& x)
   {
      switch(N)
      {
         case(0): return (0) * lambda() + u(x) * gamma();
         case(1): return (0) * lambda() + u(x) * gamma();
         case(2): return (-2) * lambda() + u(x) * gamma();
         case(3): return (-6 * x) * lambda() + u(x) * gamma();
         case(4): return (-12 * x * x) * lambda() + u(x) * gamma();
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

   double u(const double& x)
   {
      switch(N)
      {
         case(0): return 2.0;
         case(1): return x;
         case(2): return x * x;
         case(3): return x * x * x;
         case(4): return x * x * x * x;
      };
   }
};