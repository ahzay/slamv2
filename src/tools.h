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
double angle_diff(double angle1, double angle2);

template<typename Container, typename Compare, typename EqualTo>
typename Container::const_iterator
unique_min(const Container &container, Compare comp, EqualTo eq) {
    auto minIt = std::min_element(container.begin(), container.end(), comp);
    int count = std::count_if(container.begin(), container.end(),
                              [minVal = *minIt, eq](const auto &val) { return eq(val, minVal); }
    );

    if (count == 1) {
        return minIt;
    } else {
        return container.end();
    }
}

template<typename Container, typename Compare, typename EqualTo>
typename Container::const_iterator
unique_max(const Container &container, Compare comp, EqualTo eq) {
    auto maxIt = std::max_element(container.begin(), container.end(), comp);
    int count = std::count_if(container.begin(), container.end(),
                              [maxVal = *maxIt, eq](const auto &val) { return eq(val, maxVal); }
    );

    if (count == 1) {
        return maxIt;
    } else {
        return container.end();
    }
}

template<typename Container, typename Compare>
std::vector<typename Container::const_iterator>
all_max_elements(const Container &container, Compare comp) {
    auto maxIt = std::max_element(container.begin(), container.end(), comp);
    std::vector<typename Container::const_iterator> maxIts;

    if (maxIt != container.end()) {
        for (auto it = container.begin(); it != container.end(); ++it) {
            if (!comp(*it, *maxIt) && !comp(*maxIt, *it)) { // equal values
                maxIts.push_back(it);
            }
        }
    }
    return maxIts;
}


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
