//
// Created by user on 7/10/23.
//

#ifndef MAPPERNODE_TOOLS_H
#define MAPPERNODE_TOOLS_H

#include <Eigen/Core>
#include <iostream>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

using namespace Eigen;
using namespace std;

bool is_psd(const MatrixXd &A);

double vector_cov(VectorXd x, VectorXd y);

VectorXd
atan2(VectorXd y, VectorXd x);


int
sgn(double val);

enum State {
    s_begin,      // 0
    s_continuous, // 1
    // s_association,       // 2
    s_fsm,               // 3 -> 2
    s_fsm_flexible,      // 3a -> 2a ..
    s_fsm_strict,        // 3b
    s_fsm_least_squares, // 3c
    s_fsm_close,         // 3d
    s_fsm_sink,          // 3e
    s_merge              // 4 -> 3
};

#endif //MAPPERNODE_TOOLS_H
