#include <iostream>
#include "NonLinearBVP.h"

int main()
{
   NonLinearBVP nl = NonLinearBVP();
   nl.test = Test(1, 0);

   nl.ReadFormGrid("regions.txt");
   nl.ReadFormtimeGrid("time.txt");

   /*ofstream fout("result.txt");
   nl.SimpleIterations(0.0, 1e-14, 1e-14, 100, fout);
   fout.close();*/

   nl.ExplicitScheme("result.txt");
   //nl.StrangeTest();

   int asd = 1111;
}