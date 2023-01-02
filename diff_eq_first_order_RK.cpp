// Markus Buchholz, 2023
// https://en.wikipedia.org/wiki/Van_der_Pol_oscillator
// g++ diff_eq_first_order_RK.cpp -o t -I/usr/include/python3.8 -lpython3.8

#include <iostream>
#include <vector>
#include <tuple>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//--------------------------------------------------------------------------------
float MU = 5.5;
float dt = 0.01;

//--------------------------------------------------------------------------------
// dot x
float function1(float x1, float x2)
{
    return x2;
}

//--------------------------------------------------------------------------------
// dot y

float function2(float x1, float x2)
{
    return MU * (1 - x1 * x1) * x2 - x1;
}

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>> methodRungeKutta1Diff()
{

    std::vector<float> diffEq1;
    std::vector<float> diffEq2;

    std::vector<float> time;

    // init values
    float x1 = 1.5; // x1
    float x2 = 0.2; // x2
    float t = 0.0;  // init time

    diffEq1.push_back(x1);
    diffEq2.push_back(x2);
    time.push_back(t);

    for (int ii = 0; ii < 10000; ii++)
    {
        t = t + dt;
        float k11 = function1(x1, x2);
        float k12 = function2(x1, x2);

        float k21 = function1(x1 + dt / 2 * k11, x2 + dt / 2 * k12);
        float k22 = function2(x1 + dt / 2 * k11, x2 + dt / 2 * k12);

        float k31 = function1(x1 + dt / 2 * k21, x2 + dt / 2 * k22);
        float k32 = function2(x1 + dt / 2 * k21, x2 + dt / 2 * k22);

        float k41 = function1(x1 + dt * k31, x2 + dt * k32);
        float k42 = function2(x1 + dt * k31, x2 + dt * k32);

        x1 = x1 + dt / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        x2 = x2 + dt / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);

        diffEq1.push_back(x1);
        diffEq2.push_back(x2);
        time.push_back(t);
    }

    return std::make_tuple(diffEq1, diffEq2, time);
}

//---------------------------------------------------------------------------------------------------------

void plot2D(std::vector<float> xX, std::vector<float> yY)
{
    plt::title("Solution of First-Order Differential Equation ");
    plt::named_plot("solution", xX, yY);
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::legend();
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::show();
}

//---------------------------------------------------------------------------------------------------------

int main()
{
    auto vdp = methodRungeKutta1Diff();
    auto xX = std::get<0>(vdp);
    auto yY = std::get<1>(vdp);
    plot2D(xX, yY);
}