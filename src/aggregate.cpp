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
    _data_vector.back().normalize();
    if (_pose != data._pose)
        change_referential(data._pose);
}


void Aggregate::flush(const Data &data) {
    flush();
    push_back(data);

}


void Aggregate::flush() {
    //_data_vector.clear();
    // new flush works with age
    for (auto &d: _data_vector)
        d.life--;
    erase_if(_data_vector, [](auto d) { return d.life < 1; });
    //_pose.setConstant(NAN);
}

/*Aggregate::Aggregate(const MatrixX2d &mat, const Vector3d &pose) {
    _pose = pose;
    for (int i = 0; i < mat.rows(); i++) {
        _data_vector.emplace_back(mat(i, all), _pose);
    }
}*/

MatrixX2d Aggregate::get_xy_mat() const {
    MatrixX2d mat;
    mat.conservativeResize(_data_vector.size(), 2);
    for (int i = 0; i < _data_vector.size(); i++)
        mat.row(i) = _data_vector[i].get_xy();
    return mat;
}

void Aggregate::self_sort() {
    sort(_data_vector.begin(), _data_vector.end(),
         [](const auto &d1, const auto &d2) -> bool {
             return d1._measurement(1) < d2._measurement(1);
         });
}

void Aggregate::change_referential(Vector3d new_pose) {
    for (auto &d: _data_vector)
        d.change_referential(new_pose);
    _pose = new_pose;
    // sort
    self_sort();
}

void Aggregate::clear() {
    _data_vector.clear();
    _pose.setConstant(NAN);
}

void Aggregate::push_back(const Aggregate &a) {
    for (const auto &d: a._data_vector)
        push_back(d);
}

void Aggregate::reorder() {
    if(_data_vector.empty()) return;
    self_sort();
    auto max_difference_iter = _data_vector.begin();
    double max_difference = (_data_vector.front().get_xy() - _data_vector.back().get_xy()).norm();

    for (auto it = _data_vector.begin() + 1; it != _data_vector.end(); ++it) {
        double difference = (it->get_xy() - (it - 1)->get_xy()).norm();
        if (difference > max_difference) {
            max_difference = difference;
            max_difference_iter = it;
        }
    }
    rotate(_data_vector.begin(), max_difference_iter, _data_vector.end());
}

Aggregate::Aggregate() = default;
