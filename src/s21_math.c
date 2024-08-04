#include "s21_math.h"

long double s21_factorial(int n) {
  long double res;
  if (n == 0) {
    res = 1.000000;
  } else if (n < 0)
    res = S21_NAN;
  else {
    res = n * s21_factorial(n - 1);
  }
  return res;
}

long double s21_pow_int(double base, long long int exp) {
  long double result = 1;
  double sign = 1.0;
  if (exp < 0) {
    sign = -1;
    exp = -exp;
  }
  for (long long int i = 0; i < exp; i++) {
    if (sign < 0)
      result = result / base;
    else
      result = result * base;
  }
  return result;
}

int s21_abs(int x) {
  int result;
  result = x;
  if (x < 0) {
    result = -result;
  }
  return result;
}

long double s21_fabs(double n) {
  if (n < 0.0) {
    n = -n;
  }
  if (1.0 / n == S21_double_NEG_INF) n = 0.0;

  return n;
}
long double s21_fmod(double x, double y) {
  x = (long double)x;
  y = (long double)y;
  long double result = 0.0;
  int fr = 1;
  if (S21_ISNAN(x) == 1 || S21_ISNAN(y) == 1 ||
      (s21_fabs(y) == S21_double_POS_INF &&
       s21_fabs(x) == S21_double_POS_INF)) {
    result = S21_NAN;
  } else if (s21_fabs(y) == S21_double_POS_INF &&
             s21_fabs(x) != S21_double_POS_INF) {
    result = x;
  } else {
    if (x < 0) {
      x *= -1;
      fr = -1;
    }
    if (y < 0) {
      y *= -1;
    }
    result = x - s21_floor(x / y) * y;
    result *= fr;
  }
  return result;
}

long double s21_floor(double x) {
  long double result;
  if (s21_fabs(x) < s21_pow_int(2, 52) && S21_ISNAN(x) == 0 &&
      x < S21_double_POS_INF && x > S21_double_NEG_INF) {
    result = (long long int)x;
    if (s21_fabs(x - result) > 0. && s21_fabs(x) > 0.) {
      if (x < 0.) {
        result -= 1;
      }
    }
  } else {
    result = x;
  }
  return result;
}

long double s21_ceil(double x) {
  long double result;
  if (s21_fabs(x) < s21_pow_int(2, 52) && S21_ISNAN(x) == 0 &&
      x < S21_double_POS_INF && x > S21_double_NEG_INF) {
    result = (long long int)x;
    if (s21_fabs(x) > 0. && x != result) {
      if (x > 0.) {
        result += 1;
      }
    }
  } else {
    result = x;
  }
  return result;
}

long double s21_log(double n) {
  n = (long double)n;
  long double result = 0.0;
  _Bool neg_flag = 0;
  _Bool inf_flag = 0;
  _Bool zero_flag = 0;
  int power_adjust = 0;
  double add = 0.0, iteration = 1.0, series = n;
  if (n < 0.0) {
    neg_flag = 1;
  }
  if (n == 0.0) {
    zero_flag = 1;
  }
  if (n == S21_double_POS_INF) {
    inf_flag = 1;
  }
  if (zero_flag == 0 && neg_flag == 0 && inf_flag == 0) {
    while (n > 1.0) {
      n /= S21_E;
      power_adjust++;
    }
    while (n < .25) {
      n *= S21_E;
      power_adjust--;
    }
    n -= 1.0;
    series = n;
    for (int i = 1; i <= TAYLOR_ITERATIONS; i++) {
      add += series * iteration / (long double)i;
      series *= n;
      iteration = -iteration;
    }
  }
  result = add + power_adjust;
  if (zero_flag == 1) {
    result = S21_double_NEG_INF;
  }
  if (inf_flag == 1) {
    result = S21_double_POS_INF;
  }
  if (neg_flag == 1) {
    result = S21_NAN;
  }
  return result;
}

long double s21_exp(double n) {
  long double add = 1.000000;
  long double result = 1.000000;
  long double iteration = 1.000000;

  int iteration_counter = 0;
  _Bool inf_flag = 0;
  n = (long double)n;
  if (n == S21_double_NEG_INF) {
    result = 0.0;
    inf_flag = 1;
  }
  if (n < -16.0 && n != S21_double_NEG_INF) {
    result = 0.0;
    inf_flag = 1;
  }

  while (inf_flag == 0 && s21_fabs(add) > DBL_EPSILON && result < DBL_MAX &&
         result > S21_double_NEG_INF && iteration_counter < TAYLOR_ITERATIONS) {
    add *= n / iteration;
    iteration += 1.0;
    result += add;
    iteration_counter++;

    if (s21_fabs((double)result) >= DBL_MAX) {
      result = S21_double_POS_INF;
      break;
    }
  }

  return result;
}

long double s21_pow(double base, double exp) {
  long double result = 0.0;
  _Bool neg_flag = 0;
  _Bool inf_flag = 0;
  _Bool zero_flag = 0;
  _Bool nan_flag = 0;

  if (S21_ISNAN(base) == 1 || S21_ISNAN(exp) == 1) {
    nan_flag = 1;
  }
  if (base < +0.0) {
    base = s21_fabs(base);
    neg_flag = 1;
  }
  if (base == +0.0) {
    zero_flag = 1;
    if (exp < 0) {
      if (1.0 / base == S21_double_NEG_INF && exp == -1.0)
        result = S21_double_NEG_INF;
      else
        result = S21_double_POS_INF;

    } else {
      if (1.0 / base == S21_double_NEG_INF && exp == 1.0)
        result = -0.0;
      else
        result = 0.0;
    }
  }

  if (base == S21_double_POS_INF) {
    if (exp < DBL_EPSILON) {
      result = 0.0;

    } else {
      result = base;
    }
    if (((int64_t)exp) != exp) {
      neg_flag = 0;
    }
    inf_flag = 1;
  }
  if (exp == 0.0) {
    result = 1.0;
  } else {
    if (inf_flag == 0 && zero_flag == 0) {
      result = s21_exp(exp * s21_log(base));
    }
  }
  if (neg_flag == 1 && (s21_fabs(exp) >= 2.0 || s21_fabs(exp) == 1.0) &&
      ((int64_t)exp % 2 != 0)) {
    result = result * -1.0;
    if (zero_flag == 1) {
      result = -0.0;
    }
  }
  if (((int64_t)exp) != exp && neg_flag == 1 && inf_flag == 0) {
    if (zero_flag == 0 && exp != S21_double_NEG_INF) {
      result = S21_NAN;
      if (exp == S21_double_POS_INF) {
        result = S21_double_POS_INF;
      }
    } else {
      result = 0.0;
    }
  }
  if (nan_flag == 1 && exp != 0.0) {
    result = S21_NAN;
  }
  if ((base == 1.0 && neg_flag == 0) ||
      (s21_fabs(base) == 1.0 && s21_fabs(exp) == S21_double_POS_INF)) {
    result = 1.0;
  }
  return result;
}

long double s21_sqrt(double n) {
  long double result = 0.0;
  if (n > 0.0 && n > DBL_EPSILON) {
    result = s21_pow(n, 0.5);
  }
  if (S21_ISNAN(n) == 1 || n < 0) {
    result = S21_NAN;
  }
  return result;
}

long double s21_sin(double x) {
  long double res = 0.0;
  double sign = 0.0;
  x = s21_fmod(x, 2.0 * S21_PI);
  for (long long int i = 0; i < TAYLOR_ITERATIONS; i++) {
    sign = (i % 2 == 0 || i == 0) ? 1 : -1;
    res += sign * s21_pow_int(x, 2 * i + 1) / s21_factorial(2 * i + 1);
  }
  return res;
}

long double s21_cos(double x) {
  x = (long double)x;
  long double res = S21_NAN;
  if (x != S21_double_POS_INF || x != S21_double_NEG_INF || !S21_ISNAN(x))
    res = s21_sin(S21_PI_2 - x);
  return res;
}

long double s21_tan(double x) {
  x = (long double)x;
  long double res;
  long double cos = 0.0;
  if (x == S21_PI_2) {
    res = ARCTAN_PI_2;
  } else if (x == -S21_PI_2) {
    res = -(ARCTAN_PI_2);
  } else if (x == 0) {
    res = 0;
  } else {
    x = s21_fmod(x, S21_PI);
    cos = s21_cos(x);
    if (cos == 0.0) {
      res = S21_NAN;
    } else
      res = s21_sin(x) / cos;
  }
  return res;
}
long double s21_atan(double x) {
  x = (long double)x;
  long double result = 0;
  if (S21_ISNAN(x)) {
    result = S21_NAN;
  } else if (x == 1) {
    result = S21_PI_4;
  } else if (x == -1) {
    result = -S21_PI_4;
  } else if (x == S21_PI_2) {
    result = 1.003884821853887214L;
  } else if (x == -S21_PI_2) {
    result = -1.003884821853887214L;
  } else if (x == S21_double_POS_INF || x == S21_double_NEG_INF) {
    if (x < 0) {
      result = -S21_PI_2;
    } else {
      result = S21_PI_2;
    }
  } else if (-1. < x && x < 1.) {
    for (register int i = 0; i < 5000; i++) {
      result += s21_pow(-1, i) * s21_pow(x, 1 + (2 * i)) / (1 + (2 * i));
    }
  } else {
    for (register int i = 0; i < 7000; i++) {
      result += s21_pow(-1, i) * s21_pow(x, -1 - (2 * i)) / (1 + (2 * i));
    }
    result = S21_PI * s21_sqrt(x * x) / (2 * x) - result;
  }
  return result;
}
long double s21_asin(double x) {
  x = (long double)x;
  long double result = 0.0;
  if (x == 1.) {
    result = S21_PI_2;
  } else if (x == -1.) {
    result = -S21_PI_2;
  } else if (s21_fabs(x) < 1e-9) {
    result = 0.0;
  } else if (x == S21_SQRT1OVER2) {
    result = S21_PI_4;
  } else if (x == -S21_SQRT1OVER2) {
    result = -S21_PI_4;
  } else if (-1. < x && x < 1.) {
    result = s21_atan(x / s21_sqrt(1 - x * x));
  } else {
    result = S21_NAN;
  }
  return result;
}
long double s21_acos(double x) {
  x = (long double)x;
  long double result = 0.0;
  if (x == 1.) {
    result = 0.0;
  } else if (x == -1.) {
    result = S21_PI;
  } else if (x == 0.0) {
    result = S21_PI_2;
  } else if (x == S21_SQRT_05) {
    result = S21_PI_4;
  } else if (x == -S21_SQRT_05) {
    result = 3 * S21_PI_4;
  } else if (0. < x && x < 1.0) {
    result = s21_atan(s21_sqrt(1 - x * x) / x);
  } else if (-1. < x && x < 0) {
    result = S21_PI + s21_atan(s21_sqrt(1 - x * x) / x);
  } else {
    result = S21_NAN;
  }
  return result;
}
