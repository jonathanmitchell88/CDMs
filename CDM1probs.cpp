//CDM1probs.cpp
#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>

// [[Rcpp::export]]
std::vector<double> CDM1probsRcpp(const std::vector<double>& pars, const std::string& permutation) {
    double g = pars[0];
    double x1 = pars[1];
    double x2 = pars[2];
    double x3 = pars[3];
    double x4 = pars[4];
    double x5 = pars[5];

    double a0011 = x4 * x5;
    double a0101 = x2 * x3 * x4;
    double a0110 = x2 * x3 * x5;
    double a0111 = a0011 * x2 * x3;
    double a1001 = x1 * x2 * x4;
    double a1010 = x1 * x2 * x5;
    double a1100 = x1 * x3;
    double a1101 = a1001 * x3;
    double a1110 = a1010 * x3;
    double a1111 = a0011 * a1100;	
    double d = a1101 * x5;
	double a1011 = std::sqrt(a0011 * a1001 * a1010);

    double t1 = g * g;
    double t2 = 1 - t1;

    std::vector<double> q = {1, g, g, t1 + t2 * a0011, g, t1 + t2 * a0101, t1 + t2 * a0110,
                             g * (t1 + t2 * (a0011 + a0101 + a0110 - 2 * a0111)), g, t1 + t2 * a1001,
                             t1 + t2 * a1010, g * (t1 + t2 * (a0011 + a1001 + a1010 - 2 * a1011)),
                             t1 + t2 * a1100, g * (t1 + t2 * (a0101 + a1001 + a1100 - 2 * a1101)),
                             g * (t1 + t2 * (a0110 + a1010 + a1100 - 2 * a1110)),
                             t1 * (t1 + t2 * (a0011 + a0101 + a0110 + a1001 + a1010 + a1100 - 2 * (a0111 + a1011 + a1101 + a1110 - 2 * d))) + std::pow(t2, 2) * a1111};

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

    std::vector<double> p(16, 0);
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            p[i] += (1.0 / 16) * h[i][j] * q[j];
        }
    }

    if (permutation == "(23)") {
        std::vector<int> indices = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};
        std::vector<double> permuted_p(16);
        for (int i = 0; i < 16; ++i) {
            permuted_p[i] = p[indices[i]];
        }
        return permuted_p;
    } else if (permutation == "(234)") {
        std::vector<int> indices = {0, 4, 1, 5, 2, 6, 3, 7, 8, 12, 9, 13, 10, 14, 11, 15};
        std::vector<double> permuted_p(16);
        for (int i = 0; i < 16; ++i) {
            permuted_p[i] = p[indices[i]];
        }
        return permuted_p;
    }

    return p;
}

