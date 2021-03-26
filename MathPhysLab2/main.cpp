#include <iostream>
#include "NonLinearBVP.h"

int main()
{
   NonLinearBVP nl = NonLinearBVP();
   nl.test = Test(4);

   nl.ReadFormGrid("regions.txt");
   nl.FormPortrait();
   nl.FillMatrix();
   nl.AccountBound();

   nl.slae.LUDecomp();
   nl.slae.ForwardSolver();
   nl.slae.BackwardSolver();

   nl.PrintSolution("result.txt");

   int asd = 1111;
}