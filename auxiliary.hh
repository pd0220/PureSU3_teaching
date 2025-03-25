#pragma once
#pragma warning(push)
#pragma warning(disable : 4127)

// usead headers and/or libraries
#include <iostream>
#include <complex>
#include <math.h>
#include <vector>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

// set namespace
using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// imaginary unit
const complex<double> I(0., 1.);
// one over square root of 3
const double ONE_OVER_SQRT3 = 1. / sqrt(3.);
// complex 2x2 and 3x3 identity matrices
// pi
const double PI = 3.14159265359;
const Eigen::Matrix2cd ID2 = Eigen::Matrix2cd::Identity();
const Eigen::Matrix3cd ID3 = Eigen::Matrix3cd::Identity();

// Pauli matrices
const Eigen::Matrix2cd sigma1{{0, 1}, {1, 0}};
const Eigen::Matrix2cd sigma2{{0, -I}, {I, 0}};
const Eigen::Matrix2cd sigma3{{1, 0}, {0, -1}};
const vector<Eigen::Matrix2cd> Pauli = {sigma1, sigma2, sigma3};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// auxiliary lambda(s)
// general square function
auto sq = [](auto const &x)
{
    return x * x;
};

// general cubic function
auto cb = [](auto const &x)
{
    return x * x * x;
};

// general quartic function
auto quart = [](auto const &x)
{
    return sq(sq(x));
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// Eigen vector to std vector
vector<int> eigen_to_std(Eigen::VectorXi const &Eigen_v)
{
    vector<int> std_v(Eigen_v.data(), Eigen_v.data() + Eigen_v.size());
    return std_v;
}

#pragma warning(pop)