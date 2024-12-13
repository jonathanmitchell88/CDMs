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
    double a0111 = a0011 * x2 * x3;
    double a1001 = x1 * x2 * x4;
    double a1010 = x1 * x2 * x5;
    double a1100 = x1 * x3;
    double a1101 = a1001 * x3;
    double a1110 = a1010 * x3;
    double a1111 = a0011 * a1100;	
    double d = a1101 * x5;
	double a1011 = std::sqrt(a0011 * a1001 * a1010);

    double t1 = g*g;
    double t2 = 1-t1;

    std::vector<double> q = {1,g,g,t1+t2*a0011,g,t1+t2*a0101,t1+t2*a0110,g*(t1+t2*(a0011+a0101+a0110-2*a0111)),g,t1+t2*a1001,t1+t2*a1010,g*(t1+t2*(a0011+a1001+a1010-2*a1011)),t1+t2*a1100,g*(t1+t2*(a0101+a1001+a1100-2*a1101)),g*(t1+t2*(a0110+a1010+a1100-2*a1110)),t1*(t1+t2*(a0011+a0101+a0110+a1001+a1010+a1100-2*(a0111+a1011+a1101+a1110-2*d)))+t2*t2*a1111};

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

    std::vector<double> p(16);
    for(int i = 0; i < 16; ++i) {
        for(int j = 0; j < 16; ++j) {
            p[i] += h[i][j] * q[j];
        }
        p[i] /= 16;
    }

    if(permutation == "(34)") {
        std::vector<int> idx = {0,2,1,3,4,6,5,7,8,10,9,11,12,14,13,15};
        std::vector<double> temp(16);
        for(int i = 0; i < 16; ++i) temp[i] = p[idx[i]];
        p = temp;
    } else if(permutation == "(23)") {
        std::vector<int> idx = {0,1,4,5,2,3,6,7,8,9,12,13,10,11,14,15};
        std::vector<double> temp(16);
        for(int i = 0; i < 16; ++i) temp[i] = p[idx[i]];
        p = temp;
    } else if(permutation == "(243)") {
        std::vector<int> idx = {0,4,1,5,2,6,3,7,8,12,9,13,10,14,11,15};
        std::vector<double> temp(16);
        for(int i = 0; i < 16; ++i) temp[i] = p[idx[i]];
        p = temp;
    } else if(permutation == "(234)") {
        std::vector<int> idx = {0,2,4,6,1,3,5,7,8,10,12,14,9,11,13,15};
        std::vector<double> temp(16);
        for(int i = 0; i < 16; ++i) temp[i] = p[idx[i]];
        p = temp;
    } else if(permutation == "(24)") {
        std::vector<int> idx = {0,4,2,6,1,5,3,7,8,12,10,14,9,13,11,15};
        std::vector<double> temp(16);
        for(int i = 0; i < 16; ++i) temp[i] = p[idx[i]];
        p = temp;
    }

    if(*std::min_element(p.begin(), p.end()) == 0) {
        double sum_data = 0;
        for(double d : data) sum_data += d;
        for(double& pi : p) {
            if(pi == 0) pi = 1/(2*sum_data);
        }
        double sum_p = 0;
        for(double pi : p) sum_p += pi;
        for(double& pi : p) pi /= sum_p;
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
 
    double da0011dx4 = x5*t2;
    double da0011dx5 = x4*t2;

    double da0101dx2 = x3*x4*t2;
    double da0101dx3 = x2*x4*t2;
    double da0101dx4 = x2*x3*t2;

    double da0110dx2 = x3*x5*t2;
    double da0110dx3 = x2*x5*t2;
    double da0110dx5 = x2*x3*t2;

    double da1001dx1 = x2*x4*t2;
    double da1001dx2 = x1*x4*t2;
    double da1001dx4 = x1*x2*t2;

    double da1010dx1 = x2*x5*t2;
    double da1010dx2 = x1*x5*t2;
    double da1010dx5 = x1*x2*t2;

    double da1100dx1 = x3*t2;
    double da1100dx3 = x1*t2;

    double da0111dx2 = -2*a0011*x3*t2;
    double da0111dx3 = -2*a0011*x2*t2;
    double da0111dx4 = -2*x2*x3*x5*t2;
    double da0111dx5 = -2*x2*x3*x4*t2;

    double da1011dx1 = -2*a0011*x2*t2;
    double da1011dx2 = -2*a0011*x1*t2;
    double da1011dx4 = -2*a1010*t2;
    double da1011dx5 = -2*a1001*t2;

    double da1101dx1 = -2*x2*x3*x4*t2;
    double da1101dx2 = -2*a1100*x4*t2;
    double da1101dx3 = -2*a1001*t2;
    double da1101dx4 = -2*a1100*x2*t2;

    double da1110dx1 = -2*x2*x3*x5*t2;
    double da1110dx2 = -2*a1100*x5*t2;
    double da1110dx3 = -2*a1010*t2;
    double da1110dx5 = -2*a1100*x2*t2;

    double da1111dx1 = a0011*x3;
    double da1111dx3 = a0011*x1;
    double da1111dx4 = a1100*x5;
    double da1111dx5 = a1100*x4;

    double dddx2 = a0011*x1*x3;

    double t3 = 2*g;
    double t4 = 3*t1;
    double t5 = 1-t4;

    double dq0011dg = t3*(1-a0011);
    double dq0101dg = t3*(1-a0101);
    double dq0110dg = t3*(1-a0110);
    double dq0111dg = t4+t5*(a0011+a0101+a0110-2*a0111);
    double dq1001dg = t3*(1-a1001);
    double dq1010dg = t3*(1-a1010);
    double dq1011dg = t4+t5*(a0011+a1001+a1010-2*a1011);
    double dq1100dg = t3*(1-a1100);
    double dq1101dg = t4+t5*(a0101+a1001+a1100-2*a1101);
    double dq1110dg = t4+t5*(a0110+a1010+a1100-2*a1110);
    double dq1111dg = 2*(2*g*g*g+g*(1-2*g*g)*(a0011+a0101+a0110+a1001+a1010+a1100-2*(a0111+a1011+a1101+a1110-2*d))-t2*t3*a1111);

    double dq1011dx1 = g*(da1001dx1+da1010dx1+da1011dx1);
    double dq1101dx1 = g*(da1001dx1+da1100dx1+da1101dx1);
    double dq1110dx1 = g*(da1010dx1+da1100dx1+da1110dx1);
    double dq1111dx1 = t1*(da1001dx1+da1010dx1+da1100dx1+da1011dx1+da1101dx1+da1110dx1+4*t2*a0111)+t2*t2*da1111dx1;

    double dq0111dx2 = g*(da0101dx2+da0110dx2+da0111dx2);
    double dq1011dx2 = g*(da1001dx2+da1010dx2+da1011dx2);
    double dq1101dx2 = g*(da0101dx2+da1001dx2+da1101dx2);
    double dq1110dx2 = g*(da0110dx2+da1010dx2+da1110dx2);
    double dq1111dx2 = t1*(da0101dx2+da0110dx2+da1001dx2+da1010dx2+da0111dx2+da1011dx2+da1101dx2+da1110dx2+4*t2*dddx2);

    double dq0111dx3 = g*(da0101dx3+da0110dx3+da0111dx3);
    double dq1101dx3 = g*(da0101dx3+da1100dx3+da1101dx3);
    double dq1110dx3 = g*(da0110dx3+da1100dx3+da1110dx3);
    double dq1111dx3 = t1*(da0101dx3+da0110dx3+da1100dx3+da0111dx3+da1101dx3+da1110dx3+4*t2*a1011)+t2*t2*da1111dx3;

    double dq0111dx4 = g*(da0011dx4+da0101dx4+da0111dx4);
    double dq1011dx4 = g*(da0011dx4+da1001dx4+da1011dx4);
    double dq1101dx4 = g*(da0101dx4+da1001dx4+da1101dx4);
    double dq1111dx4 = t1*(da0011dx4+da0101dx4+da1001dx4+da0111dx4+da1011dx4+da1101dx4+4*t2*a1110)+t2*t2*da1111dx4;

    double dq0111dx5 = g*(da0011dx5+da0110dx5+da0111dx5);
    double dq1011dx5 = g*(da0011dx5+da1010dx5+da1011dx5);
    double dq1110dx5 = g*(da0110dx5+da1010dx5+da1110dx5);
    double dq1111dx5 = t1*(da0011dx5+da0110dx5+da1010dx5+da0111dx5+da1011dx5+da1110dx5+4*t2*a1101)+t2*t2*da1111dx5;

    std::vector<double> dpdg(17);
    std::vector<double> dpdx1(17);
    std::vector<double> dpdx2(17);
    std::vector<double> dpdx3(17);
    std::vector<double> dpdx4(17);
    std::vector<double> dpdx5(17);

    for (int i = 0; i < 16; ++i) {
        dpdg[i] = dpdq[i][1] + dpdq[i][2] + dpdq[i][3] * dq0011dg + dpdq[i][4] + dpdq[i][5] * dq0101dg + dpdq[i][6] * dq0110dg + dpdq[i][7] * dq0111dg + dpdq[i][8] + dpdq[i][9] * dq1001dg + dpdq[i][10] * dq1010dg + dpdq[i][11] * dq1011dg + dpdq[i][12] * dq1100dg + dpdq[i][13] * dq1101dg + dpdq[i][14] * dq1110dg + dpdq[i][15] * dq1111dg;
        dpdx1[i] = dpdq[i][9] * da1001dx1 + dpdq[i][10] * da1010dx1 + dpdq[i][11] * dq1011dx1 + dpdq[i][12] * da1100dx1 + dpdq[i][13] * dq1101dx1 + dpdq[i][14] * dq1110dx1 + dpdq[i][15] * dq1111dx1;
        dpdx2[i] = dpdq[i][5] * da0101dx2 + dpdq[i][6] * da0110dx2 + dpdq[i][7] * dq0111dx2 + dpdq[i][9] * da1001dx2 + dpdq[i][10] * da1010dx2 + dpdq[i][11] * dq1011dx2 + dpdq[i][13] * dq1101dx2 + dpdq[i][14] * dq1110dx2 + dpdq[i][15] * dq1111dx2;
        dpdx3[i] = dpdq[i][5] * da0101dx3 + dpdq[i][6] * da0110dx3 + dpdq[i][7] * dq0111dx3 + dpdq[i][12] * da1100dx3 + dpdq[i][13] * dq1101dx3 + dpdq[i][14] * dq1110dx3 + dpdq[i][15] * dq1111dx3;
        dpdx4[i] = dpdq[i][3] * da0011dx4 + dpdq[i][5] * da0101dx4 + dpdq[i][7] * dq0111dx4 + dpdq[i][9] * da1001dx4 + dpdq[i][11] * dq1011dx4 + dpdq[i][13] * dq1101dx4 + dpdq[i][15] * dq1111dx4;
        dpdx5[i] = dpdq[i][3] * da0011dx5 + dpdq[i][6] * da0110dx5 + dpdq[i][7] * dq0111dx5 + dpdq[i][10] * da1010dx5 + dpdq[i][11] * dq1011dx5 + dpdq[i][14] * dq1110dx5 + dpdq[i][15] * dq1111dx5;
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
