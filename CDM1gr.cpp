//CDM1gr.cpp
#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// [[Rcpp::export]]
std::vector<double> CDM1grRcpp(const std::vector<double>& pars, const std::vector<double>& data, const std::string& permutation) {
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

    double t1 = g * g;
    double t2 = 1 - t1;

    std::vector<double> q = {
        1, g, g, t1 + t2 * a0011, g, t1 + t2 * a0101, t1 + t2 * a0110,
        g * (t1 + t2 * (a0011 + a0101 + a0110 - 2 * a0111)), g, t1 + t2 * a1001,
        t1 + t2 * a1010, g * (t1 + t2 * (a0011 + a1001 + a1010 - 2 * a1011)),
        t1 + t2 * a1100, g * (t1 + t2 * (a0101 + a1001 + a1100 - 2 * a1101)),
        g * (t1 + t2 * (a0110 + a1010 + a1100 - 2 * a1110)),
        t1 * (t1 + t2 * (a0011 + a0101 + a0110 + a1001 + a1010 + a1100 - 2 * (a0111 + a1011 + a1101 + a1110 - 2 * d))) + t2 * t2 * a1111
    };

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

    std::vector<double> p = std::vector<double>(16, 0);
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            p[i] += h[i][j] * q[j];
        }
        p[i] /= 16;
    }

    if (permutation == "(23)") {
        std::vector<double> p_new = {
            p[0], p[1], p[4], p[5], p[2], p[3], p[6], p[7],
            p[8], p[9], p[12], p[13], p[10], p[11], p[14], p[15]
        };
        p = p_new;
    } else if (permutation == "(234)") {
        std::vector<double> p_new = {
            p[0], p[4], p[1], p[5], p[2], p[6], p[3], p[7],
            p[8], p[12], p[9], p[13], p[10], p[14], p[11], p[15]
        };
        p = p_new;
    }

    if (*std::min_element(p.begin(), p.end()) == 0) {
        for (auto& val : p) {
            if (val == 0) val = 1.0 / (2 * std::accumulate(data.begin(), data.end(), 0.0));
        }
        double sum_p = std::accumulate(p.begin(), p.end(), 0.0);
        for (auto& val : p) {
            val /= sum_p;
        }
    }

    std::vector<double> dlogLdp(16);
    for (int i = 0; i < 16; ++i) {
        dlogLdp[i] = -2 * data[i] / p[i];
    }

    std::vector<std::vector<double>> dpdq(16, std::vector<double>(16, 0));
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            dpdq[i][j] = h[i][j] / 16;
        }
    }

    double da0011dx4 = x5 * t2;
    double da0011dx5 = x4 * t2;

    double da0101dx2 = x3 * x4 * t2;
    double da0101dx3 = x2 * x4 * t2;
    double da0101dx4 = x2 * x3 * t2;

    double da0110dx2 = x3 * x5 * t2;
    double da0110dx3 = x2 * x5 * t2;
    double da0110dx5 = x2 * x3 * t2;

    double da1001dx1 = x2 * x4 * t2;
    double da1001dx2 = x1 * x4 * t2;
    double da1001dx4 = x1 * x2 * t2;

    double da1010dx1 = x2 * x5 * t2;
    double da1010dx2 = x1 * x5 * t2;
    double da1010dx5 = x1 * x2 * t2;

    double da1100dx1 = x3 * t2;
    double da1100dx3 = x1 * t2;

    double da0111dx2 = -2 * a0011 * x3 * t2;
    double da0111dx3 = -2 * a0011 * x2 * t2;
    double da0111dx4 = -2 * a0110 * t2;
    double da0111dx5 = -2 * a0101 * t2;

    double da1011dx1 = -2 * a0011 * x2 * t2;
    double da1011dx2 = -2 * a0011 * x1 * t2;
    double da1011dx4 = -2 * a1010 * t2;
    double da1011dx5 = -2 * a1001 * t2;

    double da1101dx1 = -2 * a0101 * t2;
    double da1101dx2 = -2 * a1100 * x4 * t2;
    double da1101dx3 = -2 * a1001 * t2;
    double da1101dx4 = -2 * a1100 * x2 * t2;

    double da1110dx1 = -2 * a0110 * t2;
    double da1110dx2 = -2 * a1100 * x5 * t2;
    double da1110dx3 = -2 * a1010 * t2;
    double da1110dx5 = -2 * a1100 * x2 * t2;

    double da1111dx1 = a0011 * x3;
    double da1111dx3 = a0011 * x1;
    double da1111dx4 = a1100 * x5;
    double da1111dx5 = a1100 * x4;

    double dddx1 = a0011 * x2 * x3;
    double dddx2 = a1111;
    double dddx3 = a0011 * x1 * x2;
    double dddx4 = a1110;
    double dddx5 = a1101;

    double t3 = 2 * g;
    double t4 = 3 * t1;
    double t5 = 1 - t4;

    double dq0011dg = t3 * (1 - a0011);
    double dq0101dg = t3 * (1 - a0101);
    double dq0110dg = t3 * (1 - a0110);
    double dq0111dg = t4 + t5 * (a0011 + a0101 + a0110 - 2 * a0111);
    double dq1001dg = t3 * (1 - a1001);
    double dq1010dg = t3 * (1 - a1010);
    double dq1011dg = t4 + t5 * (a0011 + a1001 + a1010 - 2 * a1011);
    double dq1100dg = t3 * (1 - a1100);
    double dq1101dg = t4 + t5 * (a0101 + a1001 + a1100 - 2 * a1101);
    double dq1110dg = t4 + t5 * (a0110 + a1010 + a1100 - 2 * a1110);
    double dq1111dg = 2 * (2 * std::pow(g, 3) + g * (1 - 2 * std::pow(g, 2)) * (a0011 + a0101 + a0110 + a1001 + a1010 + a1100 - 2 * (a0111 + a1011 + a1101 + a1110 - 2 * d)) - t2 * t3 * a1111);

    double dq1011dx1 = g * (da1001dx1 + da1010dx1 + da1011dx1);
    double dq1101dx1 = g * (da1001dx1 + da1100dx1 + da1101dx1);
    double dq1110dx1 = g * (da1010dx1 + da1100dx1 + da1110dx1);
    double dq1111dx1 = t1 * (da1001dx1 + da1010dx1 + da1100dx1 + da1011dx1 + da1101dx1 + da1110dx1 + 4 * t2 * dddx1) + t2 * t2 * da1111dx1;

    double dq0111dx2 = g * (da0101dx2 + da0110dx2 + da0111dx2);
    double dq1011dx2 = g * (da1001dx2 + da1010dx2 + da1011dx2);
    double dq1101dx2 = g * (da0101dx2 + da1001dx2 + da1101dx2);
    double dq1110dx2 = g * (da0110dx2 + da1010dx2 + da1110dx2);
    double dq1111dx2 = t1 * (da0101dx2 + da0110dx2 + da1001dx2 + da1010dx2 + da0111dx2 + da1011dx2 + da1101dx2 + da1110dx2 + 4 * t2 * dddx2);

    double dq0111dx3 = g * (da0101dx3 + da0110dx3 + da0111dx3);
    double dq1101dx3 = g * (da0101dx3 + da1100dx3 + da1101dx3);
    double dq1110dx3 = g * (da0110dx3 + da1100dx3 + da1110dx3);
    double dq1111dx3 = t1 * (da0101dx3 + da0110dx3 + da1100dx3 + da0111dx3 + da1101dx3 + da1110dx3 + 4 * t2 * dddx3) + t2 * t2 * da1111dx3;

    double dq0111dx4 = g * (da0011dx4 + da0101dx4 + da0111dx4);
    double dq1011dx4 = g * (da0011dx4 + da1001dx4 + da1011dx4);
    double dq1101dx4 = g * (da0101dx4 + da1001dx4 + da1101dx4);
    double dq1111dx4 = t1 * (da0011dx4 + da0101dx4 + da1001dx4 + da0111dx4 + da1011dx4 + da1101dx4 + 4 * t2 * dddx4) + t2 * t2 * da1111dx4;

    double dq0111dx5 = g * (da0011dx5 + da0110dx5 + da0111dx5);
    double dq1011dx5 = g * (da0011dx5 + da1010dx5 + da1011dx5);
    double dq1110dx5 = g * (da0110dx5 + da1010dx5 + da1110dx5);
    double dq1111dx5 = t1 * (da0011dx5 + da0110dx5 + da1010dx5 + da0111dx5 + da1011dx5 + da1110dx5 + 4 * t2 * dddx5) + t2 * t2 * da1111dx5;

    std::vector<double> dpdg(16);
    for (int i = 0; i < 16; ++i) {
        dpdg[i] = dpdq[i][1] + dpdq[i][2] + dpdq[i][3] * dq0011dg + dpdq[i][4] + dpdq[i][5] * dq0101dg +
                 dpdq[i][6] * dq0110dg + dpdq[i][7] * dq0111dg + dpdq[i][8] + dpdq[i][9] * dq1001dg +
                 dpdq[i][10] * dq1010dg + dpdq[i][11] * dq1011dg + dpdq[i][12] * dq1100dg +
                 dpdq[i][13] * dq1101dg + dpdq[i][14] * dq1110dg + dpdq[i][15] * dq1111dg;
    }

    std::vector<double> dpdx1(16);
    for (int i = 0; i < 16; ++i) {
        dpdx1[i] = dpdq[i][9] * da1001dx1 + dpdq[i][10] * da1010dx1 + dpdq[i][11] * dq1011dx1 +
                  dpdq[i][12] * da1100dx1 + dpdq[i][13] * dq1101dx1 + dpdq[i][14] * dq1110dx1 +
                  dpdq[i][15] * dq1111dx1;
    }

    std::vector<double> dpdx2(16);
    for (int i = 0; i < 16; ++i) {
        dpdx2[i] = dpdq[i][5] * da0101dx2 + dpdq[i][6] * da0110dx2 + dpdq[i][7] * dq0111dx2 +
                  dpdq[i][9] * da1001dx2 + dpdq[i][10] * da1010dx2 + dpdq[i][11] * dq1011dx2 +
                  dpdq[i][13] * dq1101dx2 + dpdq[i][14] * dq1110dx2 + dpdq[i][15] * dq1111dx2;
    }

    std::vector<double> dpdx3(16);
    for (int i = 0; i < 16; ++i) {
        dpdx3[i] = dpdq[i][5] * da0101dx3 + dpdq[i][6] * da0110dx3 + dpdq[i][7] * dq0111dx3 +
                  dpdq[i][12] * da1100dx3 + dpdq[i][13] * dq1101dx3 + dpdq[i][14] * dq1110dx3 +
                  dpdq[i][15] * dq1111dx3;
    }

    std::vector<double> dpdx4(16);
    for (int i = 0; i < 16; ++i) {
        dpdx4[i] = dpdq[i][3] * da0011dx4 + dpdq[i][5] * da0101dx4 + dpdq[i][7] * dq0111dx4 +
                  dpdq[i][9] * da1001dx4 + dpdq[i][11] * dq1011dx4 + dpdq[i][13] * dq1101dx4 +
                  dpdq[i][15] * dq1111dx4;
    }

    std::vector<double> dpdx5(16);
    for (int i = 0; i < 16; ++i) {
        dpdx5[i] = dpdq[i][3] * da0011dx5 + dpdq[i][6] * da0110dx5 + dpdq[i][7] * dq0111dx5 +
                  dpdq[i][10] * da1010dx5 + dpdq[i][11] * dq1011dx5 + dpdq[i][14] * dq1110dx5 +
                  dpdq[i][15] * dq1111dx5;
    }

    std::vector<double> derivs(6);
    derivs[0] = std::inner_product(dlogLdp.begin(), dlogLdp.end(), dpdg.begin(), 0.0) / std::accumulate(data.begin(), data.end(), 0.0);
    derivs[1] = std::inner_product(dlogLdp.begin(), dlogLdp.end(), dpdx1.begin(), 0.0) / std::accumulate(data.begin(), data.end(), 0.0);
    derivs[2] = std::inner_product(dlogLdp.begin(), dlogLdp.end(), dpdx2.begin(), 0.0) / std::accumulate(data.begin(), data.end(), 0.0);
    derivs[3] = std::inner_product(dlogLdp.begin(), dlogLdp.end(), dpdx3.begin(), 0.0) / std::accumulate(data.begin(), data.end(), 0.0);
    derivs[4] = std::inner_product(dlogLdp.begin(), dlogLdp.end(), dpdx4.begin(), 0.0) / std::accumulate(data.begin(), data.end(), 0.0);
	derivs[5] = std::inner_product(dlogLdp.begin(), dlogLdp.end(), dpdx5.begin(), 0.0) / std::accumulate(data.begin(), data.end(), 0.0);

    return derivs;
}

