#ifndef S21_MATH_H_
#define S21_MATH_H_

#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define S21_PRECISION 1.0e-6
#define S21_EPS 1.0e-7
#define S21_DBL_MAX 1.79769e+308
#define S21_double_POS_INF +(1.0 / 0.0)
#define S21_double_NEG_INF -(1.0 / 0.0)
#define S21_NAN 0.0 / 0.0
#define S21_E 2.718281828459045235
#define TAYLOR_ITERATIONS 100
#define S21_ISNAN(X) ((X) != (X))
#define S21_PI 3.14159265358979323846
#define ARCTAN_PI_2 16331239353195370.0L
#define S21_PI_2 1.570796326794896619
#define S21_SQRT_05 0.7071067811865475244
#define S21_PI_4 0.7853981633974480L
#define S21_SQRT1OVER2 0.7071067811865475244

long double s21_factorial(int n);
long double s21_fabs(double n);
long double s21_log(double n);
long double s21_exp(double n);
long double s21_pow(double base, double exp);
long double s21_sqrt(double n);
long double s21_tan(double x);
long double s21_sin(double x);
long double s21_fmod(double x, double y);
int s21_abs(int x);
long double s21_floor(double x);
long double s21_cos(double x);
long double s21_ceil(double x);
long double s21_atan(double x);
long double s21_asin(double x);
long double s21_acos(double x);

#endif  // S21_MATH_H_
