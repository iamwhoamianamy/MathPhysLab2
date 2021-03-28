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
         case(0): return (0) * lambda(u(x)) + u(x) * gamma();
         case(1): return (0) * lambda(u(x)) + u(x) * gamma();
         case(2): return (-2) * lambda(u(x)) + u(x) * gamma();
         case(3): return (-6 * x) * lambda(u(x)) + u(x) * gamma();
         case(4): return (-12 * x * x) * lambda(u(x)) + u(x) * gamma();
      };
   }

   double lambda(const double& u)
   {
      return u;
   }

   /*double lambda_elem(const vector<double>& q, const vector<double>& x_elem, const double& x)
   {
      return lambda(u_elem(q, x_elem[0])) * 2 * (x - 0.5) * (x - 1) +
             lambda(u_elem(q, x_elem[1])) * (-4 * x) * (x - 1) +
             lambda(u_elem(q, x_elem[2])) * 2 * x * (x - 0.5);
   }*/

   double u_elem(const vector<double>& q, const double& x)
   {
      return q[0] * 2 * (x - 0.5) * (x - 1) +
             q[1] * (-4 * x) * (x - 1) +
             q[2] * 2 * x * (x - 0.5);
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