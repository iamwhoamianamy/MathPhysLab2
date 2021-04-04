#include <iostream>
#include "NonLinearBVP.h"

int main()
{
   NonLinearBVP nl = NonLinearBVP();
   nl.test = Test(1, 0);

   nl.ReadFormGrid("regions.txt");
   nl.ReadFormtimeGrid("time.txt");


   nl.ExplicitScheme("result.txt");


   int asd = 1111;
}