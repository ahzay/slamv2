//
// Created by user on 7/10/23.
//

#include "aggregate.h"

MatrixX2d Aggregate::get_measurement_mat() const {
    MatrixX2d mat;
    mat.conservativeResize(_data_vector.size(), 2);
    for (int i = 0; i < _data_vector.size(); i++)
        mat.row(i) = _data_vector[i]._measurement;
    return mat;
}

MatrixX2d Aggregate::get_rotated_measurement_mat() const {
    MatrixX2d mat;
    mat.conservativeResize(_data_vector.size(), 2);
    for (int i = 0; i < _data_vector.size(); i++)
        mat.row(i) = _data_vector[i].get_rotated_measurement();
    return mat;
}

void Aggregate::push_back(Data data) {
    _data_vector.push_back(data);
}


void Aggregate::flush(Data data) {
    _data_vector.clear();
    _data_vector.push_back(data);
}


void Aggregate::flush() {
    _data_vector.clear();
}

Aggregate::Aggregate(const MatrixX2d &mat, const Vector3d &pose) {
    _pose = pose;
    for (int i = 0; i < mat.rows(); i++) {
        _data_vector.emplace_back(mat(i, all), _pose);
    }
}

MatrixX2d Aggregate::get_xy_mat() const {
    const auto dan = get_rotated_measurement_mat();
    VectorXd d = dan.col(0);
    VectorXd ori = dan.col(1);
    ArrayX2d ans;
    ans.col(0) = d.array() * ori.array().cos() + _pose(0);
    ans.col(1) = d.array() * ori.array().sin() + _pose(1);
    return ans;
}
