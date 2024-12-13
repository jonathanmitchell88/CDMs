//CDM1start1.cpp
#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

// [[Rcpp::export]]
std::vector<double> CDM1start1Rcpp(std::vector<double>& data, const std::string& permutation, double tol) {
    if (permutation == "(23)") {
        std::vector<double> temp = {data[0], data[1], data[4], data[5], data[2], data[3], data[6], data[7],
                                    data[8], data[9], data[12], data[13], data[10], data[11], data[14], data[15]};
        data = temp;
    } else if (permutation == "(234)") {
        std::vector<double> temp = {data[0], data[2], data[4], data[6], data[1], data[3], data[5], data[7],
                                    data[8], data[10], data[12], data[14], data[9], data[11], data[13], data[15]};
        data = temp;
    }

    double sum = 0;
    for (double val : data) sum += val;
    std::vector<double> x(data.size());
    for (size_t i = 0; i < data.size(); ++i) x[i] = data[i] / sum;

    auto sum_range = [&](const std::vector<int>& indices) {
        double sum = 0;
        for (int idx : indices) sum += x[idx];
        return std::max(tol, sum);
    };

    double x00oa = sum_range({0, 1, 2, 3});
    double x01oa = sum_range({4, 5, 6, 7});
    double x10oa = sum_range({8, 9, 10, 11});
    double x11oa = sum_range({12, 13, 14, 15});

    double x00ob = sum_range({0, 1, 4, 5});
    double x01ob = sum_range({2, 3, 6, 7});
    double x10ob = sum_range({8, 9, 12, 13});
    double x11ob = sum_range({10, 11, 14, 15});

    double x00oc = sum_range({0, 2, 4, 6});
    double x01oc = sum_range({1, 3, 5, 7});
    double x10oc = sum_range({8, 10, 12, 14});
    double x11oc = sum_range({9, 11, 13, 15});

    double x00ab = sum_range({0, 1, 8, 9});
	double x01ab = sum_range({2, 3, 10, 11});
	double x10ab = sum_range({4, 5, 12, 13});
    double x11ab = sum_range({6, 7, 14, 15});

    double x00bc = sum_range({0, 4, 8, 12});
    double x01bc = sum_range({1, 5, 9, 13});
    double x10bc = sum_range({2, 6, 10, 14});
    double x11bc = sum_range({3, 7, 11, 15});

    double g1 = x00oa - x11oa;
    double g2 = x00ob - x11ob;
    double g3 = x00oc - x11oc;
    double g4 = x00ab - x11ab;
    double g5 = x00bc - x11bc;

    double g = std::max(-1 + tol, std::min((g1 + g2 + g3 + g4 + g5) / 5, 1 - tol));

    auto calc_x = [&](double x00, double x11, double x01, double x10) {
        return std::max((4 * x00 * x11 - std::pow(x01 + x10, 2)) / ((2 * x00 + x01 + x10) * (x01 + x10 + 2 * x11)), tol);
    };

    double xoa = calc_x(x00oa, x11oa, x01oa, x10oa);
    double xob = calc_x(x00ob, x11ob, x01ob, x10ob);
    double xoc = calc_x(x00oc, x11oc, x01oc, x10oc);
	double xab = calc_x(x00ab, x11ab, x01ab, x10ab);
    double xbc = calc_x(x00bc, x11bc, x01bc, x10bc);

    double x3a = std::min(std::sqrt(xoa*xab/xob),1.0);
    double x5a = std::min(std::sqrt(xob*xbc/xoc),1.0);
	double x1a = std::min(xoa/x3a,1.0);
	double x2a = std::min(xab/(x3a*x5a),1.0);
	double x4a = std::min(xbc/x5a,1.0);

    auto bound = [tol](double val) { return std::max(tol, std::min(val, 1 - tol)); };

    double x1 = bound(x1a);
    double x2 = bound(x2a);
    double x3 = bound(x3a);
    double x4 = bound(x4a);
    double x5 = bound(x5a);

    return {g, x1, x2, x3, x4, x5};
}

