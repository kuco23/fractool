#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <string>
#include <fstream>
#include "../include/julia_polynomial_series.h"

using std::complex;
using std::vector;
using std::string;
using std::ofstream;
using std::endl;
using std::max;

#define M_PI 3.14159265358979323846
const complex<double> i (0, 1);

inline complex<double> horner(
  const complex_polynomial &coefs, 
  const complex<double> &z
) {
  complex<double> sum (0, 0);
  for (const complex<double> &a_i : coefs)
    sum = sum * z + a_i;
  return sum;
}

// theoretic convergence limit for a specific polynomial
double radiusJulia(const complex_polynomial &p) {
  int n = p.size();
  double an = abs(p[0]);
  double sum = 0;
  for (int i = 1; i < n; i++) sum += abs(p[i]);
  return max(
    max(1.0, 2 * sum / an), 
    pow(2 * 1.0001 / an, 1 / (double) (n - 2))
  );
}

int escapetimeMandelbrot(
  const complex<double> c, 
  const int &iterlim
) {
  int count;
  complex<double> ck = c;
  for (count = 1; count < iterlim && abs(c) <= 2; count++)
    ck = ck * ck + c;
  return (count < iterlim) ? count : 0;
}

int escapetimeJulia(
  const complex<double> z,
  const complex_polynomial &p,
  const int &iterlim
) {
  complex<double> zk = z;
  float radius = radiusJulia(p);
  int count;
  for (count = 1; count < iterlim && abs(zk) < radius; count++)
    zk = horner(p, zk);
  return (count < iterlim) ? count : 0;
}

float demMandelbrot(
  const complex<double> c,
  const int &iterlim,
  const float &radius,
  const float &overflow
) {
  complex<double> ck = c;
  complex<double> dk = 1;
  for (int i = 1; i < iterlim; i++) {
    if (max(
        max(abs(ck.imag()), abs(ck.real())),
        max(abs(dk.real()), abs(dk.imag()))) 
        > overflow) break;
    dk *= 2.0 * ck;
    dk += 1.0;
    ck *= ck;
    ck += c;
  }
  const float absck = abs(ck);
  if (absck < radius) return 0;
  const float absdk = abs(dk);
  if (absdk == 0) return 0; // rarely happens
  const float estimate = log2(absck) * absck / absdk;
  return -log2(estimate);
}

float demJulia(
  const complex<double> z,
  const complex_polynomial &p,
  const complex_polynomial &dp,
  const int &iterlim,
  const float &overflow
) {
  float radius = radiusJulia(p);

  complex<double> zk = z;
  complex<double> dk = 1;
  for (int i = 1; i < iterlim; i++) {
    if (max(
        max(abs(zk.imag()), abs(zk.real())),
        max(abs(dk.real()), abs(dk.imag()))) 
        > overflow) break;
    dk = horner(dp, zk) * dk;
    zk = horner(p, zk);
  }
  const float abszk = abs(zk);
  if (abszk < radius) return 0;
  const float absdk = abs(dk);
  if (absdk == 0) return 0; // rarely happens
  const float estimate = log2(abszk) * abszk / absdk;
  return -log2(estimate);
}

[][]float demJuliaValues(
  const complex_polynomial p,
  const complex_polynomial dp,
  const int iterlim,
  const float overflow,
  const complex<double> center,
  const float radius,
  const int px
) {
  [px][px]float values;
  complex<double> z(center.real() - radius, center.imag());
  float d = 2 * radius / px;
  complex<double> dz = d;
  for (int i = 0; i < px; i++) {
    for (int j = 0; j < px; j++) {
      values[i][j] = demJulia(z, p, dp, iterlim, overflow);
    }
    z = complex<double>(center.real() - radius, z.imag());
  }
  return values;
}