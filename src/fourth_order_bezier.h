#ifndef SPK_FOURTH_ORDER_BEZIER_H
#define SPK_FOURTH_ORDER_BEZIER_H

#include <array>
#include <tuple>

using std::array;
using std::tuple;

namespace fourth_order_bezier_polynomial {


typedef double(*PathDot)(const array<double,2>& A, const array<double,2>& B);

double cw_path_dot(const array<double,2>& A, const array<double,2>& B);
double dot(const array<double,2>& A, const array<double,2>& B);

array<double,5> 
compute_coefficients(const array<double,5>& x);

double 
compute(const array<double,5>& coefficients, double u);
void 
compute(const array<double,5>& coefficients, double u, array<double,3>& c);

void 
curve(const array<double,5>& cx, const array<double,5>& cy, double u, array<double,2>& v);

void 
curve(const array<double,5>& cx, const array<double,5>& cy, double u,
        array<double,2>& v, array<double,2>& dv);

void 
curve(const array<double,5>& cx, const array<double,5>& cy, double u, 
        array<double,2>& v, array<double,2>& dv, array<double,2>& ddv);

// Newton solver; find root; compute function defined above
double 
get_parametric_coordinate_at_max(const array<double,5>& cy, double u0, double tolerance, int max_iter);

double 
cpp_predictor_2d(const array<double,5>& cx, const array<double,5>& cy, 
             const array<double,2>& xy, double umax_curve, double du);

double 
solve(const array<double,5>& cx, double u0, double x, double tolerance, int max_iter);

tuple<double,double> 
closest_point(const array<double,5>& cx, const array<double,5>& cy, 
    const array<double,2>& xy, double umax_curve, double tolerance, int max_iter, 
    PathDot path_dot = cw_path_dot);

tuple<double, array<double,5>, array<double,5>>
scale_coordinates_for_curve_height(const array<double,5>&xt, const array<double,5>&yt, double height);

}

#endif
