//
// Created by user on 7/10/23.
//

#include "model.h"
#include "entity.h"

#include <fstream>

Model::Model(const string &file, const CmdLineOptions &options) : _options(options) {
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
    f >> _mahalanobis_strict >> _mahalanobis_flex >> _mahalanobis_aug_min >> _mahalanobis_aug >>
      _dop_sigma;
    _parameter_mins.resize(_parameter_count);
    _parameter_maxs.resize(_parameter_count);
    for (int i = 0; i < _parameter_count; i++)
        f >> _parameter_mins(i);
    for (int i = 0; i < _parameter_count; i++)
        f >> _parameter_maxs(i);
    f >> _ap_ls_forgetting_factor;
    f >> _min_npoints_mult;
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
         << _mahalanobis_aug_min << endl
         << _mahalanobis_aug << endl
         << _dop_sigma << endl
         << _parameter_mins.transpose() << endl
         << _parameter_maxs.transpose() << endl
         << _ap_ls_forgetting_factor << endl
         << _min_npoints_mult << endl;
}

void Model::dop(Entity &e, const Aggregate &a) const {
    auto df = dfsn(e, a);
    auto J = df(all, _dfs_measurement_error_indexes);
    auto H = df(all, _dfs_parameter_indexes);
    MatrixXd En = MatrixXd::Zero(_parameter_count, _parameter_count);
    for (int i = 0; i < df.rows(); i++)
        En +=
                (H(i, all).transpose() *
                 (J(i, all) * _dop_sigma * _dop_sigma
                  * J(i, all).transpose()).inverse() *
                 H(i, all));
    En = En.inverse();
    cout << "En: " << endl << En << endl;
    /*if (!is_psd(En)) {
        throw std::runtime_error("Non positive semi-definite matrix! in dop");
    }*/
    e._E = En;
}

VectorXd Model::fsn(const Entity &e, const Aggregate &a) const {
    VectorXd ans;
    ans.resize(a._data_vector.size());
    for (int i = 0; i < a._data_vector.size(); i++) {
        ans(i) = fs(e, a._data_vector[i]);
    }
    return ans;
}

MatrixXd Model::dfsn(const Entity &e, const Aggregate &a) const {
    MatrixXd ans;
    ans.resize(a._data_vector.size(), _parameter_count + 6);
    for (int i = 0; i < a._data_vector.size(); i++) {
        ans.row(i) = dfs(e, a._data_vector[i]);
    }
    return ans;
}

void Model::ap_dop(Entity &e, const Aggregate &a) const {
    auto df = dfsn(e, a);
    auto J = df(all, _dfs_measurement_error_indexes);
    auto H = df(all, _dfs_parameter_indexes);
    auto S = H * e._E * H.transpose() + J * e.m->_W_a * J.transpose();
    for (int i = 0; i < H.rows(); i++) {
        auto K = e._E * H(i, all).transpose() / S(i, i);
        auto L = e.m->I - K * H(i, all);
        e._E = L * e._E * L.transpose() + K * J(i, all) * e.m->_W_a *
                                          J(i, all).transpose() *
                                          K.transpose();
        e._E = (e._E + e._E.transpose()) / 2;
    }
    cout << "E Eigenvalues after: " << e._E.eigenvalues().transpose() << endl;
    cout << "E: " << endl << e._E.diagonal().transpose() << endl;
    if (!is_psd(e._E)) {
        throw std::runtime_error("Non positive semi-definite matrix! in ap_dop");
    }
    // system propagation
    e._E = e._E + e.m->_Q_a;// * a._data_vector.size();
    //e._E *= (1 + _ap_ls_forgetting_factor);
}

double Model::fs(const Entity &e, const Data &d) const {
    return fss(e._p, d._pose, d.get_rotated_measurement());
}
