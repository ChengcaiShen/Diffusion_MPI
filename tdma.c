//
//  function_tdma.c
//  
//
//  Created by Chengcai on 5/9/17.
//
//

#include <stdio.h>
void function_tdma(float x[], const int N, const float a[], const float b[], float c[])
{
  int i;
  
  c[0] = c[0] / b[0];
  x[0] = x[0] / b[0];
  
  for (i = 1; i < N; i++) {
    float m = 1.0f / (b[i] - a[i] * c[i - 1]);
    c[i] = c[i] * m;
    x[i] = (x[i] - a[i] * x[i - 1]) * m;
  }
  
  for (i = N - 1; i--> 0; )
    x[i] = x[i] - c[i] * x[i + 1];
}
