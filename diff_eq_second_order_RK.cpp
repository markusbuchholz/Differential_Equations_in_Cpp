// Markus Buchholz, 2023
//
// g++ diff_eq_second_order_RK.cpp -o t -I/usr/include/python3.8 -lpython3.8
//--------------------------------------------------------------------------------
/*
The most important in RK method is to understand how to code differential equation.
For differential equation:
1. Take derivative on RIGHT side
2. Program as a function1 or ... the left side
3. Run below algorithm
4.* For second order you have to compute also first derivative (here we run one variable so we have program two equations)
*/
//--------------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <tuple>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//--------------------------------------------------------------------------------
float K = 30;
float dt = 0.01;
//--------------------------------------------------------------------------------
// y dot
float function1(float y, float y_dot)
{
    return y_dot;
}

//--------------------------------------------------------------------------------
// dot y

float function2(float y, float y_dot)
{
    return -(0.9 + 0.7) * y_dot - K * y;
}

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>> methodRungeKutta2Diff()
{

    std::vector<float> diffEq1;
    std::vector<float> diffEq2;

    std::vector<float> time;

    // init values
    float x1 = 2.0; // y
    float x2 = 0.0; // y_dot
    float t = 0.0;  // init time

    diffEq1.push_back(x1);
    diffEq2.push_back(x2);
    time.push_back(t);

    for (int ii = 0; ii < 1000; ii++)
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

    plt::plot(xX, yY);
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::show();
}

//---------------------------------------------------------------------------------------------------------
void plot2D2D(std::vector<float> time, std::vector<float> xX, std::vector<float> yY)
{

    plt::title("Solution of Second-Order Differential Equation ");
    plt::named_plot("y_dot", time, yY);
    plt::named_plot("y", time, xX);
    plt::xlabel("time");
    plt::ylabel("Y");
    plt::legend();
    plt::show();
}

//---------------------------------------------------------------------------------------------------------

int main()
{
    auto diff2 = methodRungeKutta2Diff();
    auto xX = std::get<0>(diff2);
    auto yY = std::get<1>(diff2);
    auto time = std::get<2>(diff2);
    plot2D2D(time, xX, yY);
}