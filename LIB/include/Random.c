#include <iostream>
#include <cstdlib>
#include <time.h>

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int Random(int a, int b) 
{
  static bool initialized = false; 
  if (!initialized) {
    srand(time(0)); 
    initialized = true; 
  }

  if (a < b) {
    int c = b - a + 1; 
    int v = a + rand()%c; 
    return v;
  } else {
    int c = a - b + 1;
    int v = b + rand()%c;
    return v; 
  }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double Random(void) 
{
  static bool initialized = false;
  if (!initialized) {
    srand(time(0)); 
    initialized = true; 
  }

  int denom = 1048575; // 2^20 - 1.
  double v = (rand()%denom) / (denom - 1.0); 
  return v;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double Random(double a, double b)
{
  double v = a + (b - a) * Random();
  return v; 
}
//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
