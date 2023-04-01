#ifndef __tables_hh__
#define __tables_hh__

#include <array>
#include <cmath>
#include <vector>
#include <functional>
#include <eigen3/Eigen/Core>
#include "diff_eq.hh"

std::vector<double> coefs_1 = { 1 };
std::vector<double> coefs_2 = { 3./2, -1./2 };
std::vector<double> coefs_3 = { 23./12, -4./3, 5./12 };
std::vector<double> coefs_4 = { 55./24, -59./24, 37./24, -3./8 };

ButcherTable<4> ex_table_4 (
    std::array<double, 4> {   0, 1./2, 1./2,    1}, 
    std::array<double, 4> {1./6, 1./3, 1./3, 1./6},
    std::array<std::array<double, 4>, 4>  {
        std::array<double, 4> {    0,    0,    0,    0}, 
        std::array<double, 4> { 1./2,    0,    0,    0}, 
        std::array<double, 4> {    0, 1./2,    0,    0}, 
        std::array<double, 4> {    0,    0,    1,    0}
    }
);

ButcherTable<3> im_table_3 (
    std::array<double, 3> {      0.32,  0.962963,  0.962963},
    std::array<double, 3> {  0.720046,  0.720046,  0.008391},
    std::array<std::array<double, 3>, 3>  {
        std::array<double, 3> {  0.333333,         0, -0.013333}, 
        std::array<double, 3> {  0.625153,  0.333333,  0.004477}, 
        std::array<double, 3> {  9.516331, -8.886702,  0.333333} 
    }
);

ButcherTable<2> ex_table_2 (
    std::array<double, 2> { 0, 1./2},
    std::array<double, 2> { 0, 1.}, 
    std::array<std::array<double, 2>, 2>  {
        std::array<double, 2> {    0,    0}, 
        std::array<double, 2> { 1./2,    0}, 
    }
);

#endif