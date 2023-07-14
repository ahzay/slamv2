//
// Created by user on 7/10/23.
//

#include "entity.h"

Entity::Entity(Model *model, const Eigen::VectorXd &p, const Eigen::MatrixXd &E) : m(model), _p(p), _E(E) {
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

Entity::Entity(Model *model, const Aggregate &a) : m(model) {
    _t.conservativeResize(m->_post_attributes_count);
    m->ls(*this, a);
}

bool Entity::operator==(const Entity &oe) const {
    return (m == oe.m) && (_p == oe._p) && (_E == oe._E) && (_t == oe._t);
}

void Entity::print() const {
    cout << "Entity with attributes: " << endl;
    cout << "p: " << _p.transpose() << endl;
    cout << "E: " << endl << _E << endl;
}
