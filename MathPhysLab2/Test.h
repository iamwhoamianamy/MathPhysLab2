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
      return -1 * dlambdadx(x) * dudx(x) + u_prec(x) * gamma();
   }

   double lambda(const double& u)
   {
      return 3 * u * u + 1;
   }

   double dlambdadx(const double& u)
   {
      return 6 * u;
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

   double dudx(const double& x)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return 1;
         case(2): return 2 * x;
         case(3): return 3 * x * x;
         case(4): return 4 * x * x * x;
      };
   }

   double u_prec(const double& x)
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