//
// Created by user on 7/10/23.
//

#include "model.h"
#include "entity.h"
#include "tools.h"
#include <fstream>

Model::Model(const string &file) {
    int dims = 1; // for 2d
    ifstream f(file.c_str());
    if (!f.is_open())
        throw std::runtime_error("Error opening the model file!");
    f >> _model_index >> _parameter_count >> _post_attributes_count >>
      _closed_fs;
    I = MatrixXd::Identity(_parameter_count, _parameter_count);
    _dfs_parameter_indexes.conservativeResize(_parameter_count);
    for (int i = 0; i < _parameter_count; i++)
        f >> _dfs_parameter_indexes(i);
    _dfs_measurement_error_indexes.conservativeResize(dims);
    for (int i = 0; i < dims; i++)
        f >> _dfs_measurement_error_indexes(i);
    _Q_s = MatrixXd::Zero(_parameter_count, _parameter_count);
    for (int i = 0; i < _parameter_count; i++)
        f >> _Q_s(i, i);
    _Q_s = _Q_s.array().pow(2);
    _W_s = MatrixXd::Zero(dims, dims);
    for (int i = 0; i < dims; i++)
        f >> _W_s(i, i);
    _W_s = _W_s.array().pow(2);
    _Q_a = MatrixXd::Zero(_parameter_count, _parameter_count);
    for (int i = 0; i < _parameter_count; i++)
        f >> _Q_a(i, i);
    _Q_a = _Q_a.array().pow(2);
    _W_a = MatrixXd::Zero(dims, dims);
    for (int i = 0; i < dims; i++)
        f >> _W_a(i, i);
    _W_a = _W_a.array().pow(2);
    f >> _mahalanobis_strict >> _mahalanobis_flex >> _mahalanobis_aug >>
      _dop_sigma;
    for (int i = 0; i < _parameter_count; i++)
        f >> _parameter_mins(i);
    for (int i = 0; i < _parameter_count; i++)
        f >> _parameter_maxs(i);
    f.close();
    // print
    cout << "model" << endl
         << _model_index << endl
         << _parameter_count << endl
         << _post_attributes_count << endl
         << _closed_fs << endl
         << _dfs_parameter_indexes.transpose() << endl
         << _dfs_measurement_error_indexes.transpose() << endl
         << _Q_s << endl
         << _W_s << endl
         << _Q_a << endl
         << _W_a << endl
         << _mahalanobis_strict << endl
         << _mahalanobis_flex << endl
         << _mahalanobis_aug << endl
         << _dop_sigma << endl
         << _parameter_mins.transpose() << endl
         << _parameter_maxs.transpose() << endl;
}

MatrixXd Model::dop(const Entity &e, const Aggregate &a) {
    const auto data = a.get_mat();
    Matrix<double, Dynamic, 12> df = dfs(e, a);
    MatrixXd Jmes = df(all, {10});
    MatrixXd H = df(all, {0, 1, 2, 3, 4, 5});
    MatrixXd En = MatrixXd::Zero(6, 6);
    for (int i = 0; i < df.rows(); i++)
        En +=
                (H(i, all).transpose() *
                 (Jmes(i, all) * _dop_sigma * _dop_sigma * Jmes(i, all).transpose()).inverse() *
                 H(i, all));
    En = En.inverse();
    cout << "En: " << endl << En << endl;
    if (!isPsd(En)) {
        throw std::runtime_error("Non positive semi-definite matrix!");
    }
    return En;
}
