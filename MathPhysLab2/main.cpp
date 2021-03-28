#include <iostream>
#include "NonLinearBVP.h"

int main()
{
   NonLinearBVP nl = NonLinearBVP();
   nl.test = Test(1);

   nl.ReadFormGrid("regions.txt");


   ofstream fout("result.txt");

   for(int i = 0; i < 100; i++)
   {
      nl.InitMatrix();
      nl.FormPortrait();
      nl.FillMatrix();
      nl.AccountBound();

      nl.slae.LUDecomp();
      nl.slae.ForwardSolver();
      nl.slae.BackwardSolver();

      nl.q = nl.slae.b;

      nl.PrintSolution(fout);
   }

   fout.close();

   int asd = 1111;
}