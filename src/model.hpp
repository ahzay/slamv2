#ifndef MODEL_HPP
#define MODEL_HPP

#include "functions.hpp"

// general model
class Model {
public:
    // functions
    Model(VectorXd (*fs)(const VectorXd p, const VectorXd pos, const MatrixXd m),
          MatrixXd (*dfs)(const VectorXd p, const VectorXd pos, const MatrixXd m),
          MatrixXd (*dop)(const VectorXd p,
                          const VectorXd pos,
                          const MatrixXd m,
                          const double sigma),
          VectorXd (*ls)(const VectorXd loc, const MatrixXd data),
          VectorXd (*ap_ls)(const VectorXd ap_p,
                            const MatrixXd ap_dop,
                            const VectorXd loc,
                            const MatrixXd data,
                            const double sigma_d),
          bool (*safety)(VectorXd &p),
          VectorXd (*dst)(const VectorXd p, const VectorXd pos, const MatrixXd m),
          void (*init_post)(const VectorXd p, VectorXd &t, const Aggregate a),
          bool (*associate)(const VectorXd p, VectorXd &t, const Data d),
          void (*augment_post)(const VectorXd p, VectorXd &t, const Data d),
          string file) {
        _fs = fs;
        _dfs = dfs;
        _dop = dop;
        _ls = ls;
        _ap_ls = ap_ls;
        _safety = safety;
        _dst = dst;
        _init_post = init_post;
        _associate = associate;
        _augment_post = augment_post;
        /* !!file structure!!
         * model_index
         * parameter_count
         * post_attributes_count
         * closed_fs
         * dfs_parameter_indexes (vector)
         * dfs_measurement_error_indexes (vector)
         * Q_s (diagonal, sqrt)
         * W_s (diagonal, sqrt)
         * Q_a (diagonal, sqrt)
         * W_a (diagonal, sqrt)
         * mahalanobis_strict
         * mahalanobis_flex
         * mahalanobis_aug
         * dop_sigma
         */
        int dims = 1; // for 2d
        ifstream f(file.c_str());
        if (f.is_open()) {
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
                 << _dop_sigma << endl;
        } else
            throw std::runtime_error("Error opening the model file!");
    }

    // attributes defining a model
    VectorXd (*_fs)(const VectorXd p, const VectorXd pos, const MatrixXd m);

    MatrixXd (*_dfs)(const VectorXd p, const VectorXd pos, const MatrixXd m);

    MatrixXd (*_dop)(const VectorXd p,
                     const VectorXd pos,
                     const MatrixXd m,
                     const double sigma);

    VectorXd (*_ls)(const VectorXd loc, const MatrixXd data);

    VectorXd (*_ap_ls)(const VectorXd ap_p,
                       const MatrixXd ap_dop,
                       const VectorXd loc,
                       const MatrixXd data,
                       const double sigma_d);

    bool (*_safety)(VectorXd &p);

    VectorXd (*_dst)(const VectorXd p, const VectorXd pos, const MatrixXd m);

    //
    void (*_init_post)(const VectorXd p, VectorXd &t, const Aggregate a);

    bool (*_associate)(const VectorXd p, VectorXd &t, const Data d);

    void (*_augment_post)(const VectorXd p, VectorXd &t, const Data d);

    //
    VectorXi _dfs_parameter_indexes, _dfs_measurement_error_indexes;
    MatrixXd _W_a, _W_s, _Q_a, _Q_s, I;
    double _mahalanobis_strict, _mahalanobis_flex, _dop_sigma,
            _mahalanobis_aug; // thresholds
    unsigned short _parameter_count, _model_index, _post_attributes_count;
    bool _closed_fs;
};

// individual entity
class Entity {
public:
    Entity(Model *model, VectorXd _p, MatrixXd _E)
            : m(model), p(_p), E(_E) {
        t.conservativeResize(m->_post_attributes_count);
    }

    double mahalanobis(Data &data) {
        double _r = -m->_fs(p, data.p, Matrix<double, 1, Dynamic>(data.d))(0);
        //
        MatrixXd df = m->_dfs(p, data.p, Matrix<double, 1, Dynamic>(data.d));
        MatrixXd H = df(all, m->_dfs_parameter_indexes);
        MatrixXd J = df(all, m->_dfs_measurement_error_indexes);
        double S = (H * E * H.transpose() + J * m->_W_s * J.transpose())(0);
        // debug
        // if (m->_model_index != 1) {
        if (0) {
            cout << "**************fs vs dst**************" << endl;
            cout << "fs: " << m->_fs(p, data.p, Matrix<double, 1, Dynamic>(data.d))(0)
                 << endl;
            cout << "dst: "
                 << m->_dst(p, data.p, Matrix<double, 1, Dynamic>(data.d))(0) << endl;
            cout << "S: " << S << endl;
            cout << "H:" << endl << H << endl;
            cout << "J:" << endl << J << endl;
            cout << "maha fs: "
                 << m->_fs(p, data.p, Matrix<double, 1, Dynamic>(data.d))(0) *
                    m->_fs(p, data.p, Matrix<double, 1, Dynamic>(data.d))(0) / S
                 << endl;
            cout << "maha dst: "
                 << m->_dst(p, data.p, Matrix<double, 1, Dynamic>(data.d))(0) *
                    m->_dst(p, data.p, Matrix<double, 1, Dynamic>(data.d))(0) / S
                 << endl;
            cout << "dst/fs: "
                 << m->_dst(p, data.p, Matrix<double, 1, Dynamic>(data.d))(0) /
                    m->_fs(p, data.p, Matrix<double, 1, Dynamic>(data.d))(0)
                 << endl;
        }
        return _r * _r / S;
    }

    Model *m;
    VectorXd p; // parameters
    VectorXd t; // post processing attributes
    MatrixXd E; // uncertainty matrix, initialized with DOP
};

// states of finite state machine and nested state machin
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

#endif // MODEL_HPP
