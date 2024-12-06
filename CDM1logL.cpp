//CDM1logL.cpp
#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// [[Rcpp::export]]
double CDM1logLRcpp(const std::vector<double>& pars, const std::vector<double>& data, const std::string& permutation) {
    double g = pars[0];
    double x1 = pars[1];
    double x2 = pars[2];
    double x3 = pars[3];
    double x4 = pars[4];
    double x5 = pars[5];

    double a0011 = x4 * x5;
    double a0101 = x2 * x3 * x4;
    double a0110 = x2 * x3 * x5;
    double a1001 = x1 * x2 * x4;
    double a1010 = x1 * x2 * x5;
    double a1100 = x1 * x3;

    double a0111 = x5 * a0101;
    double a1011 = x5 * a1001;
    double a1101 = x1 * a0101;
    double a1110 = x1 * a0110;

    double a1111 = a0011 * a1100;
    double d = x1 * a0111;

    double t1 = std::pow(g, 2);
    double t2 = 1 - t1;

    std::vector<double> q = {1, g, g, t1 + t2 * a0011, g, t1 + t2 * a0101, t1 + t2 * a0110,
                             g * (t1 + t2 * (a0011 + a0101 + a0110 - 2 * a0111)), g, t1 + t2 * a1001,
                             t1 + t2 * a1010, g * (t1 + t2 * (a0011 + a1001 + a1010 - 2 * a1011)),
                             t1 + t2 * a1100, g * (t1 + t2 * (a0101 + a1001 + a1100 - 2 * a1101)),
                             g * (t1 + t2 * (a0110 + a1010 + a1100 - 2 * a1110)),
                             t1 * (t1 + t2 * (a0011 + a0101 + a0110 + a1001 + a1010 + a1100 - 2 * (a0111 + a1011 + a1101 + a1110 - 2 * d))) + std::pow(t2, 2) * a1111};

    std::vector<std::vector<int>> h(16, std::vector<int>(16, 1));
    std::vector<int> rows = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    std::vector<std::vector<int>> cols = {
        {2, 4, 6, 8, 10, 12, 14, 16},
        {3, 4, 7, 8, 11, 12, 15, 16},
        {2, 3, 6, 7, 10, 11, 14, 15},
        {5, 6, 7, 8, 13, 14, 15, 16},
        {2, 4, 5, 7, 10, 12, 13, 15},
        {3, 4, 5, 6, 11, 12, 13, 14},
        {2, 3, 5, 8, 10, 11, 13, 16},
        {9, 10, 11, 12, 13, 14, 15, 16},
        {2, 4, 6, 8, 9, 11, 13, 15},
        {3, 4, 7, 8, 9, 10, 13, 14},
        {2, 3, 6, 7, 9, 12, 13, 16},
        {5, 6, 7, 8, 9, 10, 11, 12},
        {2, 4, 5, 7, 9, 11, 14, 16},
        {3, 4, 5, 6, 9, 10, 15, 16},
        {2, 3, 5, 8, 9, 12, 14, 15}
    };

    for (size_t i = 0; i < rows.size(); ++i) {
        for (int col : cols[i]) {
            h[rows[i] - 1][col - 1] = -1;
        }
    }

    std::vector<double> p(16);
    for (int i = 0; i < 16; ++i) {
        p[i] = 0;
        for (int j = 0; j < 16; ++j) {
            p[i] += h[i][j] * q[j];
        }
        p[i] /= 16;
    }

    if (permutation == "(23)") {
        std::vector<int> new_order = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};
        std::vector<double> temp(16);
        for (int i = 0; i < 16; ++i) {
            temp[i] = p[new_order[i]];
        }
        p = temp;
    } else if (permutation == "(234)") {
        std::vector<int> new_order = {0, 4, 1, 5, 2, 6, 3, 7, 8, 12, 9, 13, 10, 14, 11, 15};
        std::vector<double> temp(16);
        for (int i = 0; i < 16; ++i) {
            temp[i] = p[new_order[i]];
        }
        p = temp;
    }

    if (std::any_of(p.begin(), p.end(), [](double x) { return x <= 0; })) {
        return 1e6;
    } else {
        double sum_data = std::accumulate(data.begin(), data.end(), 0.0);
        double logL = 0;
        for (size_t i = 0; i < data.size(); ++i) {
            logL += -2 * data[i] * std::log(p[i]);
        }
        return logL / sum_data;
    }
}

