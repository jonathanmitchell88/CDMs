//CDM1hess.cpp
#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>

// [[Rcpp::export]]
NumericMatrix CDM1hessRcpp(const std::vector<double>& pars, const std::vector<int>& data, const std::string& permutation) {
    int n = std::accumulate(data.begin(), data.end(), 0);
    double tol = 1.0 / (2 * n);

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

    std::vector<double> q = {1, g, g, t1 + t2 * a0011, g, t1 + t2 * a0101, t1 + t2 * a0110,
                             g * (t1 + t2 * (a0011 + a0101 + a0110 - 2 * a0111)), g, t1 + t2 * a1001,
                             t1 + t2 * a1010, g * (t1 + t2 * (a0011 + a1001 + a1010 - 2 * a1011)),
                             t1 + t2 * a1100, g * (t1 + t2 * (a0101 + a1001 + a1100 - 2 * a1101)),
                             g * (t1 + t2 * (a0110 + a1010 + a1100 - 2 * a1110)),
                             t1 * (t1 + t2 * (a0011 + a0101 + a0110 + a1001 + a1010 + a1100 - 2 * (a0111 + a1011 + a1101 + a1110 - 2 * d))) + t2 * t2 * a1111};

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
        std::vector<double> p_new = {p[0], p[1], p[4], p[5], p[2], p[3], p[6], p[7], p[8], p[9], p[12], p[13], p[10], p[11], p[14], p[15]};
        p = p_new;
    } else if (permutation == "(234)") {
        std::vector<double> p_new = {p[0], p[4], p[1], p[5], p[2], p[6], p[3], p[7], p[8], p[12], p[9], p[13], p[10], p[14], p[11], p[15]};
        p = p_new;
    }

    for (int i = 0; i < 16; ++i) {
        if (p[i] == 0) {
            p[i] = tol;
        }
    }

    double sum_p = std::accumulate(p.begin(), p.end(), 0.0);
    for (int i = 0; i < 16; ++i) {
        p[i] /= sum_p;
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

    double da0011dx4 = t2 * x5;
    double da0011dx5 = t2 * x4;

    double da0101dx2 = t2 * x3 * x4;
    double da0101dx3 = t2 * x2 * x4;
    double da0101dx4 = t2 * x2 * x3;

    double da0110dx2 = t2 * x3 * x5;
    double da0110dx3 = t2 * x2 * x5;
    double da0110dx5 = t2 * x2 * x3;

    double da1001dx1 = t2 * x2 * x4;
    double da1001dx2 = t2 * x1 * x4;
    double da1001dx4 = t2 * x1 * x2;

    double da1010dx1 = t2 * x2 * x5;
    double da1010dx2 = t2 * x1 * x5;
    double da1010dx5 = t2 * x1 * x2;

    double da1100dx1 = t2 * x3;
    double da1100dx3 = t2 * x1;

    double t3 = -2 * t2;

    double da0111dx2 = t3 * x3 * a0011;
    double da0111dx3 = t3 * x2 * a0011;
    double da0111dx4 = t3 * a0110;
    double da0111dx5 = t3 * a0101;

    double da1011dx1 = t3 * x2 * a0011;
    double da1011dx2 = t3 * x1 * a0011;
    double da1011dx4 = t3 * a1010;
    double da1011dx5 = t3 * a1001;

    double da1101dx1 = t3 * a0101;
    double da1101dx2 = t3 * x4 * a1100;
    double da1101dx3 = t3 * a1001;
    double da1101dx4 = t3 * x2 * a1100;

    double da1110dx1 = t3 * a0110;
    double da1110dx2 = t3 * x5 * a1100;
    double da1110dx3 = t3 * a1010;
    double da1110dx5 = t3 * x2 * a1100;

    double da1111dx1 = x3 * a0011;
    double da1111dx3 = x1 * a0011;
    double da1111dx4 = x5 * a1100;
    double da1111dx5 = x4 * a1100;

    double dddx1 = x2 * x3 * a0011;
    double dddx3 = x1 * x2 * a0011;

    double t4 = 2 * g;
    double t5 = 3 * t1;
    double t6 = 1 - t5;
    double t7 = 2 * t1;

    double dq0011dg = t4 * (1 - a0011);
    double dq0101dg = t4 * (1 - a0101);
    double dq0110dg = t4 * (1 - a0110);
    double dq0111dg = t5 + t6 * (a0011 + a0101 + a0110 - 2 * a0111);
    double dq1001dg = t4 * (1 - a1001);
    double dq1010dg = t4 * (1 - a1010);
    double dq1011dg = t5 + t6 * (a0011 + a1001 + a1010 - 2 * a1011);
    double dq1100dg = t4 * (1 - a1100);
    double dq1101dg = t5 + t6 * (a0101 + a1001 + a1100 - 2 * a1101);
    double dq1110dg = t5 + t6 * (a0110 + a1010 + a1100 - 2 * a1110);
    double dq1111dg = t4 * (t7 + (1 - t7) * (a0011 + a0101 + a0110 + a1001 + a1010 + a1100 - 2 * (a0111 + a1011 + a1101 + a1110 - 2 * d)) - 2 * t2 * a1111);

    double dq1011dx1 = g * (da1001dx1 + da1010dx1 + da1011dx1);
    double dq1101dx1 = g * (da1001dx1 + da1100dx1 + da1101dx1);
    double dq1110dx1 = g * (da1010dx1 + da1100dx1 + da1110dx1);
    double dq1111dx1 = t1 * (da1001dx1 + da1010dx1 + da1100dx1 + da1011dx1 + da1101dx1 + da1110dx1 + 4 * t2 * dddx1) + t2 * t2 * da1111dx1;

    double dq0111dx2 = g * (da0101dx2 + da0110dx2 + da0111dx2);
    double dq1011dx2 = g * (da1001dx2 + da1010dx2 + da1011dx2);
    double dq1101dx2 = g * (da0101dx2 + da1001dx2 + da1101dx2);
    double dq1110dx2 = g * (da0110dx2 + da1010dx2 + da1110dx2);
    double dq1111dx2 = t1 * (da0101dx2 + da0110dx2 + da1001dx2 + da1010dx2 + da0111dx2 + da1011dx2 + da1101dx2 + da1110dx2 + 4 * t2 * a1111);

    double dq0111dx3 = g * (da0101dx3 + da0110dx3 + da0111dx3);
    double dq1101dx3 = g * (da0101dx3 + da1100dx3 + da1101dx3);
    double dq1110dx3 = g * (da0110dx3 + da1100dx3 + da1110dx3);
    double dq1111dx3 = t1 * (da0101dx3 + da0110dx3 + da1100dx3 + da0111dx3 + da1101dx3 + da1110dx3 + 4 * t2 * dddx3) + t2 * t2 * da1111dx3;

    double dq0111dx4 = g * (da0011dx4 + da0101dx4 + da0111dx4);
    double dq1011dx4 = g * (da0011dx4 + da1001dx4 + da1011dx4);
    double dq1101dx4 = g * (da0101dx4 + da1001dx4 + da1101dx4);
    double dq1111dx4 = t1 * (da0011dx4 + da0101dx4 + da1001dx4 + da0111dx4 + da1011dx4 + da1101dx4 + 4 * t2 * a1110) + t2 * t2 * da1111dx4;

    double dq0111dx5 = g * (da0011dx5 + da0110dx5 + da0111dx5);
    double dq1011dx5 = g * (da0011dx5 + da1010dx5 + da1011dx5);
    double dq1110dx5 = g * (da0110dx5 + da1010dx5 + da1110dx5);
    double dq1111dx5 = t1 * (da0011dx5 + da0110dx5 + da1010dx5 + da0111dx5 + da1011dx5 + da1110dx5 + 4 * t2 * a1101) + t2 * t2 * da1111dx5;

    std::vector<double> dpdg(16, 0);
    std::vector<double> dpdx1(16, 0);
    std::vector<double> dpdx2(16, 0);
    std::vector<double> dpdx3(16, 0);
    std::vector<double> dpdx4(16, 0);
    std::vector<double> dpdx5(16, 0);

    for (int i = 0; i < 16; ++i) {
		dpdg[i] = dpdq[i][1] + dpdq[i][2] + dpdq[i][3] * dq0011dg + dpdq[i][4] + dpdq[i][5] * dq0101dg + dpdq[i][6] * dq0110dg + dpdq[i][7] * dq0111dg + dpdq[i][8] + dpdq[i][9] * dq1001dg + dpdq[i][10] * dq1010dg + dpdq[i][11] * dq1011dg + dpdq[i][12] * dq1100dg + dpdq[i][13] * dq1101dg + dpdq[i][14] * dq1110dg + dpdq[i][15] * dq1111dg;
        dpdx1[i] = dpdq[i][9] * da1001dx1 + dpdq[i][10] * da1010dx1 + dpdq[i][11] * dq1011dx1 + dpdq[i][12] * da1100dx1 + dpdq[i][13] * dq1101dx1 + dpdq[i][14] * dq1110dx1 + dpdq[i][15] * dq1111dx1;
        dpdx2[i] = dpdq[i][5] * da0101dx2 + dpdq[i][6] * da0110dx2 + dpdq[i][7] * dq0111dx2 + dpdq[i][9] * da1001dx2 + dpdq[i][10] * da1010dx2 + dpdq[i][11] * dq1011dx2 + dpdq[i][13] * dq1101dx2 + dpdq[i][14] * dq1110dx2 + dpdq[i][15] * dq1111dx2;
        dpdx3[i] = dpdq[i][5] * da0101dx3 + dpdq[i][6] * da0110dx3 + dpdq[i][7] * dq0111dx3 + dpdq[i][12] * da1100dx3 + dpdq[i][13] * dq1101dx3 + dpdq[i][14] * dq1110dx3 + dpdq[i][15] * dq1111dx3;
        dpdx4[i] = dpdq[i][3] * da0011dx4 + dpdq[i][5] * da0101dx4 + dpdq[i][7] * dq0111dx4 + dpdq[i][9] * da1001dx4 + dpdq[i][11] * dq1011dx4 + dpdq[i][13] * dq1101dx4 + dpdq[i][15] * dq1111dx4;
        dpdx5[i] = dpdq[i][3] * da0011dx5 + dpdq[i][6] * da0110dx5 + dpdq[i][7] * dq0111dx5 + dpdq[i][10] * da1010dx5 + dpdq[i][11] * dq1011dx5 + dpdq[i][14] * dq1110dx5 + dpdq[i][15] * dq1111dx5;
    }

    std::vector<double> d2logLdgdp(16);
    std::vector<double> d2logLdx1dp(16);
    std::vector<double> d2logLdx2dp(16);
    std::vector<double> d2logLdx3dp(16);
    std::vector<double> d2logLdx4dp(16);
    std::vector<double> d2logLdx5dp(16);

    for (int i = 0; i < 16; ++i) {
        d2logLdgdp[i] = 2 * data[i] / (p[i] * p[i]) * dpdg[i];
        d2logLdx1dp[i] = 2 * data[i] / (p[i] * p[i]) * dpdx1[i];
        d2logLdx2dp[i] = 2 * data[i] / (p[i] * p[i]) * dpdx2[i];
        d2logLdx3dp[i] = 2 * data[i] / (p[i] * p[i]) * dpdx3[i];
        d2logLdx4dp[i] = 2 * data[i] / (p[i] * p[i]) * dpdx4[i];
        d2logLdx5dp[i] = 2 * data[i] / (p[i] * p[i]) * dpdx5[i];
    }

    double t8 = -t4;
    double t9 = 2 * t4;
    double t10 = 2 * t4;
    double t11 = 6 * g;
    double t12 = 6 * t1;

    double d2q1111dgdd = 8 * g * (1 - 2 * t1);
    double d2q1111dgda1111 = -t2 * t9;

    double d2a0011dgdx4 = t8 * x5;
    double d2a0011dgdx5 = t8 * x4;

    double d2a0101dgdx2 = t8 * x3 * x4;
    double d2a0101dgdx3 = t8 * x2 * x4;
    double d2a0101dgdx4 = t8 * x2 * x3;

    double d2a0101dx2dx3 = t2 * x4;
    double d2a0101dx2dx4 = t2 * x3;

    double d2a0101dx3dx4 = t2 * x2;

    double d2a0110dgdx2 = t8 * x3 * x5;
    double d2a0110dgdx3 = t8 * x2 * x5;
    double d2a0110dgdx5 = t8 * x2 * x3;

    double d2a0110dx2dx3 = t2 * x5;
    double d2a0110dx2dx5 = t2 * x3;

    double d2a0110dx3dx5 = t2 * x2;

    double d2a1001dgdx1 = t8 * x2 * x4;
    double d2a1001dgdx2 = t8 * x1 * x4;
    double d2a1001dgdx4 = t8 * x1 * x2;

    double d2a1001dx1dx2 = t2 * x4;
    double d2a1001dx1dx4 = t2 * x2;

    double d2a1001dx2dx4 = t2 * x1;

    double d2a1010dgdx1 = t8 * x2 * x5;
    double d2a1010dgdx2 = t8 * x1 * x5;
    double d2a1010dgdx5 = t8 * x1 * x2;

    double d2a1010dx1dx2 = t2 * x5;
    double d2a1010dx1dx5 = t2 * x2;

    double d2a1010dx2dx5 = t2 * x1;

    double d2a1100dgdx1 = t8 * x3;
    double d2a1100dgdx3 = t8 * x1;

    double d2a0111dgdx2 = t10 * x3 * a0011;
    double d2a0111dgdx3 = t10 * x2 * a0011;
    double d2a0111dgdx4 = t10 * a0110;
    double d2a0111dgdx5 = t10 * a0101;

    double d2a0111dx2dx3 = t3 * a0011;
    double d2a0111dx2dx4 = t3 * x3 * x5;
    double d2a0111dx2dx5 = t3 * x3 * x4;

    double d2a0111dx3dx4 = t3 * x2 * x5;
    double d2a0111dx3dx5 = t3 * x2 * x4;

    double d2a0111dx4dx5 = t3 * x2 * x3;

    double d2a1011dgdx1 = t10 * x2 * a0011;
    double d2a1011dgdx2 = t10 * x1 * a0011;
    double d2a1011dgdx4 = t10 * a1010;
    double d2a1011dgdx5 = t10 * a1001;

    double d2a1011dx1dx2 = t3 * a0011;
    double d2a1011dx1dx4 = t3 * x2 * x5;
    double d2a1011dx1dx5 = t3 * x2 * x4;

    double d2a1011dx2dx4 = t3 * x1 * x5;
    double d2a1011dx2dx5 = t3 * x1 * x4;

    double d2a1011dx4dx5 = t3 * x1 * x2;

    double d2a1101dgdx1 = t10 * a0101;
    double d2a1101dgdx2 = t10 * x4 * a1100;
    double d2a1101dgdx3 = t10 * a1001;
    double d2a1101dgdx4 = t10 * x2 * a1100;

    double d2a1101dx1dx2 = t3 * x3 * x4;
    double d2a1101dx1dx3 = t3 * x2 * x4;
    double d2a1101dx1dx4 = t3 * x2 * x3;

    double d2a1101dx2dx3 = t3 * x1 * x4;
    double d2a1101dx2dx4 = t3 * a1100;

    double d2a1101dx3dx4 = t3 * x1 * x2;

    double d2a1110dgdx1 = t10 * a0110;
    double d2a1110dgdx2 = t10 * x5 * a1100;
    double d2a1110dgdx3 = t10 * a1010;
    double d2a1110dgdx5 = t10 * x2 * a1100;

    double d2a1110dx1dx2 = t3 * x3 * x5;
    double d2a1110dx1dx3 = t3 * x2 * x5;
    double d2a1110dx1dx5 = t3 * x2 * x3;

    double d2a1110dx2dx3 = t3 * x1 * x5;
    double d2a1110dx2dx5 = t3 * a1100;

    double d2a1110dx3dx5 = t3 * x1 * x2;

    double d2a1111dx1dx3 = a0011;
    double d2a1111dx1dx4 = x3 * x5;
    double d2a1111dx1dx5 = x3 * x4;

    double d2a1111dx3dx4 = x1 * x5;
    double d2a1111dx3dx5 = x1 * x4;

    double d2a1111dx4dx5 = a1100;

    double d2ddx1dx2 = x3 * a0011;
    double d2ddx1dx3 = x2 * a0011;
    double d2ddx1dx4 = a0110;
    double d2ddx1dx5 = a0101;

    double d2ddx2dx3 = x1 * a0011;
    double d2ddx2dx4 = x5 * a1100;
    double d2ddx2dx5 = x4 * a1100;

    double d2ddx3dx4 = a1010;
    double d2ddx3dx5 = a1001;

    double d2ddx4dx5 = x2 * a1100;

    double d2q0011dg2 = 2 * (1 - a0011);
    double d2q0101dg2 = 2 * (1 - a0101);
    double d2q0110dg2 = 2 * (1 - a0110);
    double d2q0111dg2 = t11 * (1 - a0011 - a0101 - a0110 + 2 * a0111);
    double d2q1001dg2 = 2 * (1 - a1001);
    double d2q1010dg2 = 2 * (1 - a1010);
    double d2q1011dg2 = t11 * (1 - a0011 - a1001 - a1010 + 2 * a1011);
    double d2q1100dg2 = 2 * (1 - a1100);
    double d2q1101dg2 = t11 * (1 - a0101 - a1001 - a1100 + 2 * a1101);
    double d2q1110dg2 = t11 * (1 - a0110 - a1010 - a1100 + 2 * a1110);
    double d2q1111dg2 = 2 * (t12 + (1 - t12) * (a0011 + a0101 + a0110 + a1001 + a1010 + a1100 - 2 * (a0111 + a1011 + a1101 + a1110 - 2 * d)) - 2 * (1 - 3 * t1) * a1111);

    double d2q1011dgdx1 = da1001dx1 + da1010dx1 + da1011dx1 + g * (d2a1001dgdx1 + d2a1010dgdx1 + d2a1011dgdx1);
    double d2q1101dgdx1 = da1001dx1 + da1100dx1 + da1101dx1 + g * (d2a1001dgdx1 + d2a1100dgdx1 + d2a1101dgdx1);
    double d2q1110dgdx1 = da1010dx1 + da1100dx1 + da1110dx1 + g * (d2a1010dgdx1 + d2a1100dgdx1 + d2a1110dgdx1);
    double d2q1111dgdx1 = t4 * da1001dx1 + t4 * da1010dx1 + t4 * da1100dx1 + t4 * da1011dx1 + t4 * da1101dx1 + t4 * da1110dx1 + d2q1111dgdd * dddx1 + d2q1111dgda1111 * da1111dx1 + t1 * (d2a1001dgdx1 + d2a1010dgdx1 + d2a1100dgdx1 + d2a1011dgdx1 + d2a1101dgdx1 + d2a1110dgdx1);

    double d2q0111dgdx2 = da0101dx2 + da0110dx2 + da0111dx2 + g * (d2a0101dgdx2 + d2a0110dgdx2 + d2a0111dgdx2);
    double d2q1011dgdx2 = da1001dx2 + da1010dx2 + da1011dx2 + g * (d2a1001dgdx2 + d2a1010dgdx2 + d2a1011dgdx2);
    double d2q1101dgdx2 = da0101dx2 + da1001dx2 + da1101dx2 + g * (d2a0101dgdx2 + d2a1001dgdx2 + d2a1101dgdx2);
    double d2q1110dgdx2 = da0110dx2 + da1010dx2 + da1110dx2 + g * (d2a0110dgdx2 + d2a1010dgdx2 + d2a1110dgdx2);
    double d2q1111dgdx2 = t4 * da0101dx2 + t4 * da0110dx2 + t4 * da1001dx2 + t4 * da1010dx2 + t4 * da0111dx2 + t4 * da1011dx2 + t4 * da1101dx2 + t4 * da1110dx2 + d2q1111dgdd * a1111 + t1 * (d2a0101dgdx2 + d2a0110dgdx2 + d2a1001dgdx2 + d2a1010dgdx2 + d2a0111dgdx2 + d2a1011dgdx2 + d2a1101dgdx2 + d2a1110dgdx2);

    double d2q0111dgdx3 = da0101dx3 + da0110dx3 + da0111dx3 + g * (d2a0101dgdx3 + d2a0110dgdx3 + d2a0111dgdx3);
    double d2q1101dgdx3 = da0101dx3 + da1100dx3 + da1101dx3 + g * (d2a0101dgdx3 + d2a1100dgdx3 + d2a1101dgdx3);
    double d2q1110dgdx3 = da0110dx3 + da1100dx3 + da1110dx3 + g * (d2a0110dgdx3 + d2a1100dgdx3 + d2a1110dgdx3);
    double d2q1111dgdx3 = t4 * da0101dx3 + t4 * da0110dx3 + t4 * da1100dx3 + t4 * da0111dx3 + t4 * da1101dx3 + t4 * da1110dx3 + d2q1111dgdd * dddx3 + d2q1111dgda1111 * da1111dx3 + t1 * (d2a0101dgdx3 + d2a0110dgdx3 + d2a1100dgdx3 + d2a0111dgdx3 + d2a1101dgdx3 + d2a1110dgdx3);

    double d2q0111dgdx4 = da0011dx4 + da0101dx4 + da0111dx4 + g * (d2a0011dgdx4 + d2a0101dgdx4 + d2a0111dgdx4);
    double d2q1011dgdx4 = da0011dx4 + da1001dx4 + da1011dx4 + g * (d2a0011dgdx4 + d2a1001dgdx4 + d2a1011dgdx4);
    double d2q1101dgdx4 = da0101dx4 + da1001dx4 + da1101dx4 + g * (d2a0101dgdx4 + d2a1001dgdx4 + d2a1101dgdx4);
    double d2q1111dgdx4 = t4 * da0011dx4 + t4 * da0101dx4 + t4 * da1001dx4 + t4 * da0111dx4 + t4 * da1011dx4 + t4 * da1101dx4 + d2q1111dgdd * a1110 + d2q1111dgda1111 * da1111dx4 + t1 * (d2a0011dgdx4 + d2a0101dgdx4 + d2a1001dgdx4 + d2a0111dgdx4 + d2a1011dgdx4 + d2a1101dgdx4);

    double d2q0111dgdx5 = da0011dx5 + da0110dx5 + da0111dx5 + g * (d2a0011dgdx5 + d2a0110dgdx5 + d2a0111dgdx5);
    double d2q1011dgdx5 = da0011dx5 + da1010dx5 + da1011dx5 + g * (d2a0011dgdx5 + d2a1010dgdx5 + d2a1011dgdx5);
    double d2q1110dgdx5 = da0110dx5 + da1010dx5 + da1110dx5 + g * (d2a0110dgdx5 + d2a1010dgdx5 + d2a1110dgdx5);
    double d2q1111dgdx5 = t4 * da0011dx5 + t4 * da0110dx5 + t4 * da1010dx5 + t4 * da0111dx5 + t4 * da1011dx5 + t4 * da1110dx5 + d2q1111dgdd * a1101 + d2q1111dgda1111 * da1111dx5 + t1 * (d2a0011dgdx5 + d2a0110dgdx5 + d2a1010dgdx5 + d2a0111dgdx5 + d2a1011dgdx5 + d2a1110dgdx5);

    double d2q1011dx1dx2 = g * (d2a1001dx1dx2 + d2a1010dx1dx2 + d2a1011dx1dx2);
    double d2q1101dx1dx2 = g * (d2a1001dx1dx2 + d2a1101dx1dx2);
    double d2q1110dx1dx2 = g * (d2a1010dx1dx2 + d2a1110dx1dx2);
    double d2q1111dx1dx2 = t1 * (d2a1001dx1dx2 + d2a1010dx1dx2 + d2a1011dx1dx2 + d2a1101dx1dx2 + d2a1110dx1dx2 + 4 * t2 * d2ddx1dx2);

    double d2q1101dx1dx3 = g * (t2 + d2a1101dx1dx3);
    double d2q1110dx1dx3 = g * (t2 + d2a1110dx1dx3);
    double d2q1111dx1dx3 = t1 * (t2 + d2a1101dx1dx3 + d2a1110dx1dx3 + 4 * t2 * d2ddx1dx3) + t2 * t2 * d2a1111dx1dx3;

    double d2q1011dx1dx4 = g * (d2a1001dx1dx4 + d2a1011dx1dx4);
    double d2q1101dx1dx4 = g * (d2a1001dx1dx4 + d2a1101dx1dx4);
    double d2q1111dx1dx4 = t1 * (d2a1001dx1dx4 + d2a1011dx1dx4 + d2a1101dx1dx4 + 4 * t2 * d2ddx1dx4) + t2 * t2 * d2a1111dx1dx4;

    double d2q1011dx1dx5 = g * (d2a1010dx1dx5 + d2a1011dx1dx5);
    double d2q1110dx1dx5 = g * (d2a1010dx1dx5 + d2a1110dx1dx5);
    double d2q1111dx1dx5 = t1 * (d2a1010dx1dx5 + d2a1011dx1dx5 + d2a1110dx1dx5 + 4 * t2 * d2ddx1dx5) + t2 * t2 * d2a1111dx1dx5;

    double d2q0111dx2dx3 = g * (d2a0101dx2dx3 + d2a0110dx2dx3 + d2a0111dx2dx3);
    double d2q1101dx2dx3 = g * (d2a0101dx2dx3 + d2a1101dx2dx3);
    double d2q1110dx2dx3 = g * (d2a0110dx2dx3 + d2a1110dx2dx3);
    double d2q1111dx2dx3 = t1 * (d2a0101dx2dx3 + d2a0110dx2dx3 + d2a0111dx2dx3 + d2a1101dx2dx3 + d2a1110dx2dx3 + 4 * t2 * d2ddx2dx3);

    double d2q0111dx2dx4 = g * (d2a0101dx2dx4 + d2a0111dx2dx4);
    double d2q1011dx2dx4 = g * (d2a1001dx2dx4 + d2a1011dx2dx4);
    double d2q1101dx2dx4 = g * (d2a0101dx2dx4 + d2a1001dx2dx4 + d2a1101dx2dx4);
    double d2q1111dx2dx4 = t1 * (d2a0101dx2dx4 + d2a1001dx2dx4 + d2a0111dx2dx4 + d2a1011dx2dx4 + d2a1101dx2dx4 + 4 * t2 * d2ddx2dx4);

    double d2q0111dx2dx5 = g * (d2a0110dx2dx5 + d2a0111dx2dx5);
    double d2q1011dx2dx5 = g * (d2a1010dx2dx5 + d2a1011dx2dx5);
    double d2q1110dx2dx5 = g * (d2a0110dx2dx5 + d2a1010dx2dx5 + d2a1110dx2dx5);
    double d2q1111dx2dx5 = t1 * (d2a0110dx2dx5 + d2a1010dx2dx5 + d2a0111dx2dx5 + d2a1011dx2dx5 + d2a1110dx2dx5 + 4 * t2 * d2ddx2dx5);

    double d2q0111dx3dx4 = g * (d2a0101dx3dx4 + d2a0111dx3dx4);
    double d2q1101dx3dx4 = g * (d2a0101dx3dx4 + d2a1101dx3dx4);
    double d2q1111dx3dx4 = t1 * (d2a0101dx3dx4 + d2a0111dx3dx4 + d2a1101dx3dx4 + 4 * t2 * d2ddx3dx4) + t2 * t2 * d2a1111dx3dx4;

    double d2q0111dx3dx5 = g * (d2a0110dx3dx5 + d2a0111dx3dx5);
    double d2q1110dx3dx5 = g * (d2a0110dx3dx5 + d2a1110dx3dx5);
    double d2q1111dx3dx5 = t1 * (d2a0110dx3dx5 + d2a0111dx3dx5 + d2a1110dx3dx5 + 4 * t2 * d2ddx3dx5) + t2 * t2 * d2a1111dx3dx5;

    double d2q0111dx4dx5 = g * (t2 + d2a0111dx4dx5);
    double d2q1011dx4dx5 = g * (t2 + d2a1011dx4dx5);
    double d2q1111dx4dx5 = t1 * (t2 + d2a0111dx4dx5 + d2a1011dx4dx5 + 4 * t2 * d2ddx4dx5) + t2 * t2 * d2a1111dx4dx5;

    std::vector<std::vector<double>> d2pdgx(16, std::vector<double>(16, 0));

    for (int i = 0; i < 16; ++i) {
        d2pdgx[0][i] = dpdq[i][3] * d2q0011dg2 + dpdq[i][5] * d2q0101dg2 + dpdq[i][6] * d2q0110dg2 + dpdq[i][7] * d2q0111dg2 + dpdq[i][9] * d2q1001dg2 + dpdq[i][10] * d2q1010dg2 + dpdq[i][11] * d2q1011dg2 + dpdq[i][12] * d2q1100dg2 + dpdq[i][13] * d2q1101dg2 + dpdq[i][14] * d2q1110dg2 + dpdq[i][15] * d2q1111dg2;
        d2pdgx[1][i] = dpdq[i][9] * d2a1001dgdx1 + dpdq[i][10] * d2a1010dgdx1 + dpdq[i][11] * d2q1011dgdx1 + dpdq[i][12] * d2a1100dgdx1 + dpdq[i][13] * d2q1101dgdx1 + dpdq[i][14] * d2q1110dgdx1 + dpdq[i][15] * d2q1111dgdx1;
        d2pdgx[2][i] = dpdq[i][5] * d2a0101dgdx2 + dpdq[i][6] * d2a0110dgdx2 + dpdq[i][7] * d2q0111dgdx2 + dpdq[i][9] * d2a1001dgdx2 + dpdq[i][10] * d2a1010dgdx2 + dpdq[i][11] * d2q1011dgdx2 + dpdq[i][13] * d2q1101dgdx2 + dpdq[i][14] * d2q1110dgdx2 + dpdq[i][15] * d2q1111dgdx2;
        d2pdgx[3][i] = dpdq[i][5] * d2a0101dgdx3 + dpdq[i][6] * d2a0110dgdx3 + dpdq[i][7] * d2q0111dgdx3 + dpdq[i][12] * d2a1100dgdx3 + dpdq[i][13] * d2q1101dgdx3 + dpdq[i][14] * d2q1110dgdx3 + dpdq[i][15] * d2q1111dgdx3;
        d2pdgx[4][i] = dpdq[i][3] * d2a0011dgdx4 + dpdq[i][5] * d2a0101dgdx4 + dpdq[i][7] * d2q0111dgdx4 + dpdq[i][9] * d2a1001dgdx4 + dpdq[i][11] * d2q1011dgdx4 + dpdq[i][13] * d2q1101dgdx4 + dpdq[i][15] * d2q1111dgdx4;
        d2pdgx[5][i] = dpdq[i][3] * d2a0011dgdx5 + dpdq[i][6] * d2a0110dgdx5 + dpdq[i][7] * d2q0111dgdx5 + dpdq[i][10] * d2a1010dgdx5 + dpdq[i][11] * d2q1011dgdx5 + dpdq[i][14] * d2q1110dgdx5 + dpdq[i][15] * d2q1111dgdx5;
        d2pdgx[6][i] = dpdq[i][9] * d2a1001dx1dx2 + dpdq[i][10] * d2a1010dx1dx2 + dpdq[i][11] * d2q1011dx1dx2 + dpdq[i][13] * d2q1101dx1dx2 + dpdq[i][14] * d2q1110dx1dx2 + dpdq[i][15] * d2q1111dx1dx2;
        d2pdgx[7][i] = dpdq[i][12] * t2 + dpdq[i][13] * d2q1101dx1dx3 + dpdq[i][14] * d2q1110dx1dx3 + dpdq[i][15] * d2q1111dx1dx3;
        d2pdgx[8][i] = dpdq[i][9] * d2a1001dx1dx4 + dpdq[i][11] * d2q1011dx1dx4 + dpdq[i][13] * d2q1101dx1dx4 + dpdq[i][15] * d2q1111dx1dx4;
        d2pdgx[9][i] = dpdq[i][10] * d2a1010dx1dx5 + dpdq[i][11] * d2q1011dx1dx5 + dpdq[i][14] * d2q1110dx1dx5 + dpdq[i][15] * d2q1111dx1dx5;
        d2pdgx[10][i] = dpdq[i][5] * d2a0101dx2dx3 + dpdq[i][6] * d2a0110dx2dx3 + dpdq[i][7] * d2q0111dx2dx3 + dpdq[i][13] * d2q1101dx2dx3 + dpdq[i][14] * d2q1110dx2dx3 + dpdq[i][15] * d2q1111dx2dx3;
        d2pdgx[11][i] = dpdq[i][5] * d2a0101dx2dx4 + dpdq[i][7] * d2q0111dx2dx4 + dpdq[i][9] * d2a1001dx2dx4 + dpdq[i][11] * d2q1011dx2dx4 + dpdq[i][13] * d2q1101dx2dx4 + dpdq[i][15] * d2q1111dx2dx4;
        d2pdgx[12][i] = dpdq[i][6] * d2a0110dx2dx5 + dpdq[i][7] * d2q0111dx2dx5 + dpdq[i][10] * d2a1010dx2dx5 + dpdq[i][11] * d2q1011dx2dx5 + dpdq[i][14] * d2q1110dx2dx5 + dpdq[i][15] * d2q1111dx2dx5;
        d2pdgx[13][i] = dpdq[i][5] * d2a0101dx3dx4 + dpdq[i][7] * d2q0111dx3dx4 + dpdq[i][13] * d2q1101dx3dx4 + dpdq[i][15] * d2q1111dx3dx4;
        d2pdgx[14][i] = dpdq[i][6] * d2a0110dx3dx5 + dpdq[i][7] * d2q0111dx3dx5 + dpdq[i][14] * d2q1110dx3dx5 + dpdq[i][15] * d2q1111dx3dx5;
        d2pdgx[15][i] = dpdq[i][3] * t2 + dpdq[i][7] * d2q0111dx4dx5 + dpdq[i][11] * d2q1011dx4dx5 + dpdq[i][15] * d2q1111dx4dx5;
    }

    std::vector<std::vector<double>> hess(6, std::vector<double>(6, 0));

    hess[0][0] = std::inner_product(d2logLdgdp.begin(), d2logLdgdp.end(), dpdg.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[0].begin(), 0.0);
    hess[1][0] = std::inner_product(d2logLdgdp.begin(), d2logLdgdp.end(), dpdx1.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[1].begin(), 0.0);
    hess[1][1] = std::inner_product(d2logLdx1dp.begin(), d2logLdx1dp.end(), dpdx1.begin(), 0.0);
    hess[2][0] = std::inner_product(d2logLdgdp.begin(), d2logLdgdp.end(), dpdx2.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[2].begin(), 0.0);
    hess[2][1] = std::inner_product(d2logLdx1dp.begin(), d2logLdx1dp.end(), dpdx2.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[6].begin(), 0.0);
    hess[2][2] = std::inner_product(d2logLdx2dp.begin(), d2logLdx2dp.end(), dpdx2.begin(), 0.0);
    hess[3][0] = std::inner_product(d2logLdgdp.begin(), d2logLdgdp.end(), dpdx3.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[3].begin(), 0.0);
    hess[3][1] = std::inner_product(d2logLdx1dp.begin(), d2logLdx1dp.end(), dpdx3.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[7].begin(), 0.0);
    hess[3][2] = std::inner_product(d2logLdx2dp.begin(), d2logLdx2dp.end(), dpdx3.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[10].begin(), 0.0);
    hess[3][3] = std::inner_product(d2logLdx3dp.begin(), d2logLdx3dp.end(), dpdx3.begin(), 0.0);
    hess[4][0] = std::inner_product(d2logLdgdp.begin(), d2logLdgdp.end(), dpdx4.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[4].begin(), 0.0);
    hess[4][1] = std::inner_product(d2logLdx1dp.begin(), d2logLdx1dp.end(), dpdx4.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[8].begin(), 0.0);
    hess[4][2] = std::inner_product(d2logLdx2dp.begin(), d2logLdx2dp.end(), dpdx4.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[11].begin(), 0.0);
    hess[4][3] = std::inner_product(d2logLdx3dp.begin(), d2logLdx3dp.end(), dpdx4.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[13].begin(), 0.0);
    hess[4][4] = std::inner_product(d2logLdx4dp.begin(), d2logLdx4dp.end(), dpdx4.begin(), 0.0);
    hess[5][0] = std::inner_product(d2logLdgdp.begin(), d2logLdgdp.end(), dpdx5.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[5].begin(), 0.0);
    hess[5][1] = std::inner_product(d2logLdx1dp.begin(), d2logLdx1dp.end(), dpdx5.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[9].begin(), 0.0);
    hess[5][2] = std::inner_product(d2logLdx2dp.begin(), d2logLdx2dp.end(), dpdx5.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[12].begin(), 0.0);
    hess[5][3] = std::inner_product(d2logLdx3dp.begin(), d2logLdx3dp.end(), dpdx5.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[14].begin(), 0.0);
    hess[5][4] = std::inner_product(d2logLdx4dp.begin(), d2logLdx4dp.end(), dpdx5.begin(), 0.0) + std::inner_product(dlogLdp.begin(), dlogLdp.end(), d2pdgx[15].begin(), 0.0);
    hess[5][5] = std::inner_product(d2logLdx5dp.begin(), d2logLdx5dp.end(), dpdx5.begin(), 0.0);

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            hess[i][j] /= n;
        }
    }
	
    NumericVector a0;
    NumericVector a1;
	NumericVector a2;
	NumericVector a3;
	NumericVector a4;
	NumericVector a5;
    NumericMatrix b;
	a0 = hess[0];
    a1 = hess[1];
	a2 = hess[2];
	a3 = hess[3];
	a4 = hess[4];
	a5 = hess[5];
    b = Rcpp::cbind(a0, a1);
	b = Rcpp::cbind(b, a2);
	b = Rcpp::cbind(b, a3);
	b = Rcpp::cbind(b, a4);
	b = Rcpp::cbind(b, a5);
    b = Rcpp::transpose(b);
    
	return b;
}

