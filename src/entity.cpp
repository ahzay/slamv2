//
// Created by user on 7/10/23.
//

#include "entity.h"

Entity::Entity(Model *model, Eigen::VectorXd p, Eigen::MatrixXd E) : m(model), _p(p), _E(E) {
    _t.conservativeResize(m->_post_attributes_count);
}

double Entity::mahalanobis(Data &data) {
    double _r = -m->fs(*this, data);
    //
    //MatrixXd df = m->_dfs(p, data.p, Matrix<double, 1, Dynamic>(data.d));
    auto df = m->dfs(*this, data);
    MatrixXd H = df(all, m->_dfs_parameter_indexes);
    MatrixXd J = df(all, m->_dfs_measurement_error_indexes);
    double S = (H * _E * H.transpose() + J * m->_W_s * J.transpose())(0);
    return _r * _r / S;
}