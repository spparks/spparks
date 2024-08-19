#ifndef SPK_VARIABLE_HAZ_H
#define SPK_VARIABLE_HAZ_H
#include <cmath>
#include <cstdio>
#include <string>

namespace variable_haz {

class VariableHaz {
public:
  /*
   * u0: parametric coordinate at maximum pool width
   * h0: user specified haz at maxium pool width
   * ht: user specified haz at tail apex
   * n: user specified shape parameter; should be n=3/2 or n=2
   */
  VariableHaz(double u0, double h0, double ht, double n)
      : u0(u0), h0(h0), ht(ht), n(n), h1((ht - h0) / (std::pow(u0, n))) {}

  void print() const {
    std::string su0 =
        "u0 = %6.4f: parametric coordinate at maximum pool width.\n";
    std::string sh0 = "h0 = %6.4f: user specified haz at maxium pool width.\n";
    std::string sht = "ht = %6.4f: user specified haz at tail apex.\n";
    std::string sn =
        "n = %3.1f: user specified shape parameter; should be n=3/2 or n=2.\n";
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
  double u0, h0, ht, n;
  double h1;
};

} // namespace variable_haz

#endif
