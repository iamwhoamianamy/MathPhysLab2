#include <iostream>
#include "NonLinearBVP.h"

int main()
{
   NonLinearBVP nl = NonLinearBVP();

   nl.ReadFormGrid("regions.txt");
   nl.FormPortrait();

   int asd = 1111;
}