//CDM1start2.cpp
#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

// [[Rcpp::export]]
std::vector<double> CDM1start2Rcpp(std::vector<double> data, const std::string& permutation, double tol) {
    if (permutation == "(23)") {
        std::vector<int> indices = {0,1,4,5,2,3,6,7,8,9,12,13,10,11,14,15};
        std::vector<double> temp(16);
        for (int i = 0; i < 16; ++i) {
            temp[i] = data[indices[i]];
        }
        data = temp;
    } else if (permutation == "(234)") {
        std::vector<int> indices = {0,2,4,6,1,3,5,7,8,10,12,14,9,11,13,15};
        std::vector<double> temp(16);
        for (int i = 0; i < 16; ++i) {
            temp[i] = data[indices[i]];
        }
        data = temp;
    }

    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    std::vector<double> x(16);
    for (int i = 0; i < 16; ++i) {
        x[i] = data[i] / sum;
    }

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

    std::vector<double> y(16, 0);
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            y[i] += h[i][j] * x[j];
        }
    }

    double g = std::max(-1+tol, std::min((y[1]+y[2]+y[4]+y[8])/4, 1-tol));

    double t1 = g * g;
    double t2 = 1 - t1;

    auto clamp = [tol](double val) { return std::max(tol, std::min(val, 1-tol)); };

    double r0011 = clamp((y[3]-t1)/t2);
    double r0101 = clamp((y[5]-t1)/t2);
    double r0110 = clamp((y[6]-t1)/t2);
    double r1001 = clamp((y[9]-t1)/t2);
    double r1100 = clamp((y[12]-t1)/t2);

    double x5 = clamp(std::sqrt(r0110*r0011/r0101));
    double x4 = clamp(r0011/x5);
    double x3 = clamp(std::sqrt(r1100*r0101/r1001));
    double x2 = clamp(r0101/(x3*x4));
    double x1 = clamp(r1100/x3);

    return {g, x1, x2, x3, x4, x5};
}

