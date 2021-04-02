#include <iostream>
#include "NonLinearBVP.h"

int main()
{
   NonLinearBVP nl = NonLinearBVP();
   nl.test = Test(1, 0);

   nl.ReadFormGrid("regions.txt");

   nl.SimpleIterations(1e-14, 1e-14, 100, "result.txt");

   int asd = 1111;
}