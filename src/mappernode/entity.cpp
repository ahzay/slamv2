//
// Created by user on 7/10/23.
//

#include "entity.h"

Entity::Entity(Model *model, const Eigen::VectorXd &p, const Eigen::MatrixXd &E) : m(model), _p(p), _E(E) {
    _t.conservativeResize(m->_post_attributes_count);
    _a = new Aggregate();
}


double Entity::mahalanobis(const Data &data) const {
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
    _a = new Aggregate(a);
    m->ls(*this, a, 0);
}

bool Entity::operator==(const Entity &oe) const {
    return (m == oe.m) && (_p == oe._p) && (_E == oe._E) && (_t == oe._t);
}

void Entity::print() const {
    cout << "Entity with attributes: " << endl;
    cout << "p: " << _p.transpose() << endl;
    cout << "E: " << endl << _E << endl;
}

Entity::~Entity() {
    delete _a;
}

void Entity::update() {
    if (_a && augmented_last_scan) {
        augmented_last_scan = false;
        auto e_old = *this;
        m->ls(*this, *_a, 1);
        _E += m->_Q_a;
        if ((e_old._p - _p).array().abs().sum() < 0.01) // didn't change
            updateLeft--;
        if (e_old.mahalanobis(*_a).first < mahalanobis(*_a).first) { // it became worse
            *this = e_old;
            updateLeft--;
        }
        // prune all points that have a maha higher than avg+std dev
        /*const auto thresh = mahalanobis(*_a).first + mahalanobis(*_a).second;
        erase_if(_a->_data_vector,
                 [*this, thresh](auto const &d) { return this->mahalanobis(d) > thresh; });*/
        //
        if (!updateLeft) {
            delete _a;
            _a = nullptr;
        }
    }
}

Entity::Entity(const Entity &other) : m(other.m), _p(other._p), _E(other._E) {
    _t = other._t;
    if (other._a) {
        _a = new Aggregate(*other._a);  // Create a deep copy of the _a object
    } else {
        _a = nullptr;
    }
}

// Define the copy assignment operator to also handle proper copying
Entity &Entity::operator=(const Entity &other) {
    if (this == &other) {
        return *this; // Self-assignment, no need to do anything
    }

    m = other.m;
    _p = other._p;
    _E = other._E;
    _t = other._t;

    // Clean up existing _a (if any)
    delete _a;

    if (other._a) {
        _a = new Aggregate(*other._a);  // Create a deep copy of the _a object
    } else {
        _a = nullptr;
    }

    return *this;
}

pair<double, double> Entity::mahalanobis(const Aggregate &aggregate) const {
    double sum_distance = 0.0;
    double sum_squared_distance = 0.0;
    for (const auto &d: aggregate._data_vector) {
        double distance = mahalanobis(d);
        sum_distance += distance;
        sum_squared_distance += distance * distance;
    }

    double average_distance = sum_distance / aggregate._data_vector.size();
    double variance = (sum_squared_distance / aggregate._data_vector.size()) - (average_distance * average_distance);
    double standard_deviation = std::sqrt(variance);

    return std::make_pair(average_distance, standard_deviation);
}


