//CDM4start.cpp
#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

// [[Rcpp::export]]
std::vector<double> CDM4startRcpp(const std::vector<double>& data, const std::string& permutation, double tol) {
    double n = 0;
    for (double d : data) n += d;

    std::vector<double> p = data;
    if (permutation == "(34)") {
        p = {data[0], data[2], data[1], data[3], data[4], data[6], data[5], data[7],
             data[8], data[10], data[9], data[11], data[12], data[14], data[13], data[15]};
    } else if (permutation == "(23)") {
        p = {data[0], data[1], data[4], data[5], data[2], data[3], data[6], data[7],
             data[8], data[9], data[12], data[13], data[10], data[11], data[14], data[15]};
    } else if (permutation == "(243)") {
        p = {data[0], data[4], data[1], data[5], data[2], data[6], data[3], data[7],
             data[8], data[12], data[9], data[13], data[10], data[14], data[11], data[15]};
    } else if (permutation == "(234)") {
        p = {data[0], data[2], data[4], data[6], data[1], data[3], data[5], data[7],
             data[8], data[10], data[12], data[14], data[9], data[11], data[13], data[15]};
    } else if (permutation == "(24)") {
        p = {data[0], data[4], data[2], data[6], data[1], data[5], data[3], data[7],
             data[8], data[12], data[10], data[14], data[9], data[13], data[11], data[15]};
    }

    for (double& d : p) d /= n;

    std::array<std::array<double, 16>, 16> h = {};
    for (int i = 0; i < 16; ++i) {
        h[i].fill(1);
    }

    h[1] = {1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1};
    h[2] = {1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1};
    h[3] = {1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1};
    h[4] = {1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1};
    h[5] = {1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1};
    h[6] = {1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1};
    h[7] = {1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1};
    h[8] = {1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1};
    h[9] = {1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1};
    h[10] = {1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1};
    h[11] = {1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1};
    h[12] = {1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1};
    h[13] = {1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1};
    h[14] = {1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1};
    h[15] = {1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1};

    std::vector<double> q(16);
    for (int i = 0; i < 16; ++i) {
        q[i] = 0;
        for (int j = 0; j < 16; ++j) {
            q[i] += h[i][j] * p[j];
        }
    }

    double g = std::max(-1 + tol, std::min((q[1] + q[2] + q[4] + q[8]) / 4, 1 - tol));

    auto clamp = [tol](double x) { return std::max(tol, std::min(x, 1 - tol)); };

    double r0011 = clamp((q[3] - g * g) / (1 - g * g));
    double r0101 = clamp((q[5] - g * g) / (1 - g * g));
    double r0110 = clamp((q[6] - g * g) / (1 - g * g));
    double r1001 = clamp((q[9] - g * g) / (1 - g * g));
    double r1010 = clamp((q[10] - g * g) / (1 - g * g));
    double r1100 = clamp((q[12] - g * g) / (1 - g * g));

    if (std::abs(g) < tol) g = g > 0 ? tol : -tol;

    double r0111 = clamp(0.5 * (r0011 + r0101 + r0110 - (q[7] - g * g * g) / (g * (1 - g * g))));
    double r1011 = clamp(0.5 * (r0011 + r1001 + r1010 - (q[11] - g * g * g) / (g * (1 - g * g))));
    double r1101 = clamp(0.5 * (r0101 + r1001 + r1100 - (q[13] - g * g * g) / (g * (1 - g * g))));
    double r1110 = clamp(0.5 * (r0110 + r1010 + r1100 - (q[14] - g * g * g) / (g * (1 - g * g))));

    double a1 = 0.5 * (r0111 + r1011 + r1101 + r1110 - 0.5 * (r0011 + r0101 + r0110 + r1001 + r1010 + r1100 - (q[15] - g * g * g * g) / (g * g * (1 - g * g))));
    double a2 = -(1 - g * g) / (4 * g * g);

    double b1 = -r0111 * r1011 * r1011 * r1101 * r1101 * r1110;
    double b2 = r1011 * r1101 * r1110 * (r0101 * r1011 - r0111 * r1001);
    double b3 = r0110 * r1001 * r1011 * r1101;
    double b4 = r0011 * r1001 * r1100;
    double b5 = r0011 * r1001 * r1101 * r1110;

    double alpha1 = a2 * b1;
    double alpha2 = a2 * b2 + a1 * b5;
    double alpha3 = a2 * b3 - b5;
    double alpha4 = a2 * b4;

    std::complex<double> c1 = std::complex<double>(-0.5, std::sqrt(3.0) / 2);
    std::complex<double> Delta0 = alpha3 * alpha3 - 3 * alpha2 * alpha4;
    std::complex<double> Delta1 = 2 * alpha3 * alpha3 * alpha3 - 9 * alpha2 * alpha3 * alpha4 + 27 * alpha1 * alpha4 * alpha4;
    std::complex<double> C = std::pow((Delta1 + std::sqrt(Delta1 * Delta1 - 4.0 * Delta0 * Delta0 * Delta0)) / 2.0, 1.0 / 3.0);

    std::vector<double> ds = {
        std::real(-1.0 / (3 * alpha4) * (alpha3 + C + Delta0 / C)),
        std::real(-1.0 / (3 * alpha4) * (alpha3 + c1 * C + Delta0 / (c1 * C))),
        std::real(-1.0 / (3 * alpha4) * (alpha3 + c1 * c1 * C + Delta0 / (c1 * c1 * C)))
    };

    double d = *std::find_if(ds.begin(), ds.end(), [](double x) { return x > 0 && x < 1; });

    double y1 = clamp(d/r0111);
    double y3 = clamp(d/std::sqrt(r0011*r1001*r1010));
    double y5 = clamp(d/r1101);
    double y6 = clamp(r1101/d*std::sqrt(r0011*r1010/r1001));
    double y8 = clamp(r0111*r1001*r1100/(r1101*d));
    double y2 = clamp(y3*r0111*r1001*r1010/(d*d));
    double y4 = clamp(y6*d*d/(r0111*r1010*r1100));
    double y7 = clamp((r0110*y4-r0111)/(y4*y8*(1-y6)));

    return {g, y1, y2, y3, y4, y5, y6, y7, y8};
}

