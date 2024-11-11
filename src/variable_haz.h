#ifndef SPK_VARIABLE_HAZ_H
#define SPK_VARIABLE_HAZ_H
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace variable_haz {

static double validate_h0_param(double h0) {
  if (h0 <= 0.0) {
    throw std::invalid_argument("Variable heat affected zone parameters "
                                "invalid.\n\t Requires h0>0.0\n");
  }
  return h0;
}

static double validate_ht_param(double h0, double ht) {
  if (h0 <= 0.0 || ht < h0) {
    throw std::invalid_argument("Variable heat affected zone parameters "
                                "invalid.\n\t Requires h0>0.0 && ht>=h0\n");
  }
  return ht;
}

static double validate_n_param(double n) {
  if (n < 1.5 || n > 2.0) {
    throw std::invalid_argument("Variable heat affected zone parameters "
                                "invalid.\n\t Requires n>1.5 && n<2.0\n");
  }
  return n;
}

static double validate_u0_param(double u0) {
  if (u0 <= 0.0 || u0 >= 1.0) {
    throw std::invalid_argument("Variable heat affected zone parameters "
                                "invalid.\n\t Requires u0>0.0 && v0<1.0\n");
  }
  return u0;
}

class VariableHaz {
public:
  /*
   * u0: parametric coordinate at maximum pool width
   * h0: user specified haz at maxium pool width
   * ht: user specified haz at tail apex
   * n: user specified shape parameter; should be n=3/2 or n=2
   */
  VariableHaz(double u0, double h0, double ht, double n)
      : u0(validate_u0_param(u0)), h0(validate_h0_param(h0)),
        ht(validate_ht_param(h0, ht)), n(validate_n_param(n)),
        h1((ht - h0) / (std::pow(u0, n))) {}

  void print() const {
    std::string su0 =
        "u0 = %6.4f: parametric coordinate at maximum pool width.\n";
    std::string sh0 = "h0 = %6.4f: user specified haz at maxium pool width.\n";
    std::string sht = "ht = %6.4f: user specified haz at tail apex.\n";
    std::string sn = "n = %3.1f: user specified shape parameter; should be "
                     "n=3/2 or n=2.\n";
    std::printf(su0.c_str(), u0);
    std::printf(sh0.c_str(), h0);
    std::printf(sht.c_str(), ht);
    std::printf(sn.c_str(), n);
  }

  double compute(double u) const {
    // Variable haz computation
    double haz = h0;
    auto du = (u0 - u);
    if (du >= 0)
      haz = h0 + h1 * std::pow(du, n);
    return haz;
  }

  double compute_mobility(double distance_from_meltpool, double u) const {
    double d = distance_from_meltpool;
    double mz;
    if (d <= 0) {
      mz = 1.0;
    } else {
      double haz = compute(u);
      if (d < haz)
        mz = 1.0 - d / haz;
      else
        mz = 0.0;
    }
    return mz;
  }

  double get_ht() const { return ht; }
  double get_h0() const { return h0; }

  double normalize_distance(double distance_from_meltpool, double u) const {
    double d = distance_from_meltpool;
    double s;
    if (d <= 0) {
      s = 0.0;
    } else {
      double haz = compute(u);
      if (d < haz)
        s = d / haz;
      else
        s = 1.0;
    }
    return s;
  }

private:
  const double u0, h0, ht, n;
  const double h1;
};

} // namespace variable_haz

#endif
