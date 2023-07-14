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

void Aggregate::push_back(const Data &data) {
    _data_vector.push_back(data);
    _pose = data._pose;
}


void Aggregate::flush(const Data &data) {
    _data_vector.clear();
    push_back(data);
}


void Aggregate::flush() {
    _data_vector.clear();
    _pose.setConstant(NAN);
}

Aggregate::Aggregate(const MatrixX2d &mat, const Vector3d &pose) {
    _pose = pose;
    for (int i = 0; i < mat.rows(); i++) {
        _data_vector.emplace_back(mat(i, all), _pose);
    }
}

MatrixX2d Aggregate::get_xy_mat() const {
    MatrixX2d mat;
    mat.conservativeResize(_data_vector.size(), 2);
    for (int i = 0; i < _data_vector.size(); i++)
        mat.row(i) = _data_vector[i].get_xy();
    return mat;
}

Aggregate::Aggregate() = default;
