//
// Created by user on 7/12/23.
//

#ifndef MAPPERNODE_EKF_H
#define MAPPERNODE_EKF_H

#include "entity.h"
#include "data.h"
class Ekf {
public:
    Ekf(const Entity &e);

    bool update(const Data &d, bool is_strict);

    double mahalanobis(const Aggregate &a) const;

    Entity _e;
    // attributes
    MatrixXd S;
    Matrix<double, 1, -1> df, H, J;
    MatrixXd K;
    double _mahalanobis{}, r{}, tol{};
};


#endif //MAPPERNODE_EKF_H
