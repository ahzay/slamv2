//
// Created by user on 7/12/23.
//

#include "ekf.h"

Ekf::Ekf(const Entity &e) : _e(e) {}

bool Ekf::update(const Data &d, const bool is_strict) {
    cout << "  in ekf" << endl;
    cout << "   E:" << endl << _e._E << endl;
    cout << "   p:" << _e._p.transpose() << endl;
    if (_e.m->safety(_e))
        cout << "forced safe values on p" << endl;
    r = -_e.m->fs(_e, d);
    df = _e.m->dfs(_e, d);
    H = df(_e.m->_dfs_parameter_indexes);
    J = df(_e.m->_dfs_measurement_error_indexes);
    cout << "   H:" << H << endl;
    cout << "   J:" << J << endl;
    S = (H * _e._E * H.transpose() + J * _e.m->_W_s * J.transpose());
    cout << "   S: " << S << endl;
    K = _e._E * H.transpose() * S.inverse();
    cout << "   K: " << K.transpose() << endl;
    _mahalanobis = (r * r * S.inverse()).diagonal().sum();
    //_mahalanobis = (dst * dst * S.inverse()).diagonal().sum();
    cout << "   mahalanobis: " << _mahalanobis << endl;
    if (is_strict)
        tol = _e.m->_mahalanobis_strict;
    else tol = _e.m->_mahalanobis_flex;
    if (_mahalanobis > tol)
        return true;
    // update
    cout << "   p" << endl;
    _e._p = _e._p + K * r;
    cout << "   E" << endl;
    _e._E = (_e.m->I - K * H) * _e._E;
    // system propagation
    _e._E = _e._E + _e.m->_Q_s;
    return false; // no bp or aberrant measure
}
