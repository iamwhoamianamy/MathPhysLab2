#include <iostream>
#include "NonLinearBVP.h"

int main()
{
   NonLinearBVP nl = NonLinearBVP();
   nl.test = Test(1, 1);

   nl.ReadFormGrid("regions.txt");
   nl.ReadFormtimeGrid("time.txt");

   // ДЛЯ НЕЛИНЕЙНОСТИ
   //ofstream fout("result.txt");
   //nl.SimpleIterations(0.0, 1e-14, 1e-14, 100, fout);
   //fout.close();

   // ДЛЯ НЕСТАЦИОНАРНОСТИ
   nl.ExplicitScheme(1e-14, 1e-14, 100, "result.txt");

   int asd = 1111;
}