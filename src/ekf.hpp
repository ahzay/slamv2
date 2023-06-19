#ifndef EKF_HPP
#define EKF_HPP
#include "model.hpp"
class Ekf { // have to update to use Data class
public:
  Ekf(const Model *m, VectorXd p, MatrixXd E) : _model(m), _p(p), _E(E) {}
  int update(const VectorXd data, const VectorXd pose, const double tol) {
    cout << "  in ekf" << endl;
    cout << "   E:" << endl << _E.diagonal().transpose() << endl;
    cout << "   p:" << _p.transpose() << endl;
    if (_model->_safety(_p))
      cout << "forced safe values on p" << endl;
    _m = data;
    // cout << "   r" << endl;
    _r = -_model->_fs(_p, pose, _m);
    // testing
    double dst = -_model->_dst(_p, pose, _m)(0);
    cout << "**************fs vs dst**************" << endl;
    cout << "fs: " << _r << endl;
    cout << "dst: " << dst << endl;
    cout << "*************************************" << endl;
    //_r(0) = dst; // comment if not using
    // cout << "   df" << endl;
    MatrixXd df = _model->_dfs(_p, pose, _m);
    H = df(all, _model->_dfs_parameter_indexes);
    cout << "   H:" << H << endl;
    J = df(all, _model->_dfs_measurement_error_indexes); // see dfsn_0
    cout << "   J:" << J << endl;
    S = (H * _E * H.transpose() + J * _model->_W_s * J.transpose());
    cout << "   S: " << S << endl;
    K = _E * H.transpose() * S.inverse();
    cout << "   K: " << K.transpose() << endl;
    //  m_dist
    // cout << "test maha: " << endl << _r * _r.transpose() * S.inverse() <<
    // endl;
    _mahalanobis = (_r * _r.transpose() * S.inverse()).diagonal().sum();
    //_mahalanobis = (dst * dst * S.inverse()).diagonal().sum();
    cout << "   mahalanobis: " << _mahalanobis << endl;
    if (_mahalanobis > tol)
      return 1;
    // update
    cout << "   p" << endl;
    _p = _p + K * _r;
    cout << "   E" << endl;
    _E = (_model->I - K * H) * _E;
    // system propagation
    _E = _E + _model->_Q_s;
    return 0; // no bp or aberrant measure
  }

  // attributes
  VectorXd _p;
  MatrixXd _E, S;
  Matrix<double, 1, Dynamic> _m;
  MatrixXd H, J, K;
  double _mahalanobis;
  VectorXd _r;
  const Model *_model;
};

class Iekf {
public:
  // E is Sigma
  Iekf(Entity *e) {
    _e = e;
    _m0.conservativeResize(4, 0); // i data size, j data num
    //_p0.conservativeResize(e->m->_parameter_count, 0); // j, data num
  }
  void add(const VectorXd data, const VectorXd pose) {
    _pose = pose;
    _m0.conservativeResize(NoChange, _m0.cols() + 1);
    _m0(all, last) = data;
    // have to pad out p0 (for each measurement I guess?), check pk
    _p0 = _e->p;
    // init for algo
    _mk = _m0;
    _pk = _p0;
  }
  int update(const double tol) {
    // ofstream o1("data1.dat");
    //  ofstream o2("data2.dat");
    const int m = _m0.cols();
    cout << "E: " << endl << _e->E.diagonal().transpose() << endl;
    /*cout << "system propagation" << endl;
    _e->E = _e->E + _e->m->_Q_a;
    cout << "E: " << endl << _e->E << endl;
    const int m = _m0.cols();
    cout << "m: " << m << endl;
    for (j = 0; j < 1; j++) {
      cout << "!!!!!!!!!!!!!!!! IEKF ITERATION!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      // jacobians
      _df = _e->m->_dfs(_pk, _pose, _mk.transpose());
      cout << "df(0): " << _df(0, all) << endl;
      H = _df(all, _e->m->_dfs_parameter_indexes);
      J = _df(all, _e->m->_dfs_measurement_error_indexes); // see dfsn_0
      _r = _e->m->_fs(_pk, _pose, _mk.transpose());
      if (_r.array().abs().mean() < 1)
        o1 << j << " " << _r.array().abs().mean() << endl;
      cout << "r: " << _r(seqN(0, min(10, m))).transpose() << endl;
      _dr.resizeLike(_r);
      for (i = 0; i < m; i++)
        _dr(i) = (J(i, all) * (_m0({ 2, 3 }, i) - _mk({ 2, 3 }, i)) +
                  H(i, all) * (_p0 - _pk))(0);

      //_dr = (J * (_m0({ 2, 3 }, all) - _mk({ 2, 3 }, all)).diagonal() +
      //       H * (_p0 - _pk));
      cout << "dr: " << _dr(seqN(0, min(10, m))).transpose() << endl;
      S = (H * _e->E * H.transpose() + J * _e->m->_W_a * J.transpose());
      // alt S
      // S.resizeLike(_dr);
      // for (i = 0; i < m; i++)
      //  S(i) = (H(i, all) * _e->E * H(i, all).transpose() +
      //          J(i, all) * _e->m->_W_a * J(i, all).transpose())(0);
      // cout << "S: " << endl << S(seqN(0, 10), seqN(0, 10)) << endl;

      // alt K
      // K = _e->E * H.transpose();
      // for (i = 0; i < m; i++)
      //  K(all, i) /= S(i);
      // cout << "K: " << endl << K(all, seqN(0, min(10, m))) << endl;
      if ((_e->E * H.transpose() * S.inverse()).array().isNaN().sum()) {
        if (j == 0) {
          return 0;
        } else {
          cout << "BREAK IEKF K NAN" << endl;
          break;
        }
      }
      K = _e->E * H.transpose() * S.inverse();
      cout << "max H: " << H.maxCoeff() << " min H: " << H.minCoeff()
           << " avg H: " << H.array().abs().mean() << endl;
      // cout << "max invH: " << H.inverse().maxCoeff()
      //      << " min invH: " << H.inverse().minCoeff()
      //      << " avg invH: " << H.inverse().array().abs().mean() << endl;
      // o2 << j << " " << (-H.inverse() * _r).array().abs().mean() << endl;
      // cout << "-invH*r: " << (-H.inverse() * _r).transpose() << endl;
      cout << "max S: " << S.maxCoeff() << " min S: " << S.minCoeff()
           << " avg S: " << S.array().abs().mean() << endl;
      cout << "max K: " << K.maxCoeff() << " min K: " << K.minCoeff()
           << " avg K: " << K.array().abs().mean() << endl;
      cout << "max r: " << _r.maxCoeff() << " min r: " << _r.minCoeff()
           << " avg r:" << _r.array().abs().mean() << endl;
      cout << "max dr: " << _dr.maxCoeff() << " min dr: " << _dr.minCoeff()
           << " avg dr:" << _dr.array().abs().mean() << endl;
      cout << "- K * (_r + _dr): " << (-K * (_r + _dr)).transpose() << endl;
      cout << "param diffs: " << (_pk - (_p0 - K * (_r + _dr)) * 1).transpose()
           << endl;
      _pk = _p0 - K * (_r + _dr) * 1;
      cout << "pk: " << _pk.transpose() << endl;
      if (_e->m->_safety(_pk)) {
        cout << "BREAK IEKF SAFETY" << endl;
        break; // force safe range for parameters, break if unsafe
      }
      Matrix<double, 2, Dynamic> G = _e->m->_W_a * J.transpose() * S.inverse();
      //Matrix<double, 2, Dynamic> G = _e->m->_W_a * J.transpose();
      //for (i = 0; i < m; i++) {
      //  G(all, i) /= S(i);
      //}
    Matrix<double, 2, Dynamic> Gg;
    Gg.conservativeResize(NoChange, _m0.cols());
    for (i = 0; i < m; i++) {
      Gg(all, i) = G(all, i) * (_r(i) + _dr(i));
    }
    // Matrix<double, Dynamic, 2> G =
    //   (_e->m->_W_a * J.transpose() *
    //    S.completeOrthogonalDecomposition().pseudoInverse())
    //     .transpose();
    // VectorXd rr = (_r + _dr);
    // auto Gg = G.array().colwise() * rr.array();
    // cout << "Gg: " << endl << Gg(all, seqN(0, min(10, m))) << endl;
    // cout << "mk: " << endl << _mk(all, seqN(0, min(10, m))) << endl;
    //_mk({ 2, 3 }, all) = _m0({ 2, 3 }, all) - Gg;
    //_mk({ 2, 3 }, all) = _m0({ 2, 3 }, all) - G * (_r + _dr);
    //  cout << "mk: " << mk.transpose() << endl;
    //  cout << "pk: " << pk.transpose() << endl;
    cout << "sum: " << _r.array().abs().sum() / _r.cols() << endl << endl;
    if (_r.array().abs().sum() / _r.cols() < tol)
      break;
  }
  cout << "j: " << j << endl;
  // jacs
  _df = _e->m->_dfs(_pk, _pose, _mk.transpose());
  H = _df(all, _e->m->_dfs_parameter_indexes);
  J = _df(all, _e->m->_dfs_measurement_error_indexes);
  //
  _e->E = _e->m->I - K * H;
  // L = _e->m->I - K * H;
  //_e->E = L * _e->E * L.transpose() +
  //         K * J * _e->m->_W_a * J.transpose() * K.transpose();
  */
    VectorXd buf =
        _e->m->_ap_ls(_e->p, _e->E, _pose, _m0.transpose(), _e->m->_dop_sigma);
    _e->p = buf(seq(0, 5));
    // cout << "mk before" << _mk.row(2)(seq(0, 10)) << endl;
    _mk.row(2) = buf(seq(6, last));
    // cout << "mk after" << _mk.row(2)(seq(0, 10)) << endl;
    _df = _e->m->_dfs(_e->p, _pose, _mk.transpose());
    H = _df(all, _e->m->_dfs_parameter_indexes);
    // cout << "H: " << endl << H << endl;
    cout << "m: " << m << endl;
    cout << "H row-rank: " << rowRank(H) << endl;
    cout << "H col-rank: " << rowRank(H.transpose()) << endl;
    J = _df(all, _e->m->_dfs_measurement_error_indexes);
    cout << "J row-rank: " << rowRank(J) << endl;
    cout << "J col-rank: " << rowRank(J.transpose()) << endl;
    S = H * _e->E * H.transpose() + J * _e->m->_W_a * J.transpose();
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(S);
    cout << "S invertible: " << qr.isInvertible() << endl;
    cout << "S determinant: " << qr.absDeterminant() << endl;
    cout << "S row-rank: " << rowRank(S) << endl;
    cout << "S col-rank: " << rowRank(S.transpose()) << endl;
    //   alt S
    // S.resize(m);
    // for (i = 0; i < m; i++)
    //  S(i) = (H(i, all) * _e->E * H(i, all).transpose() +
    //          J(i, all) * _e->m->_W_a * J(i, all).transpose())(0);
    // cout << "S cond: "
    //      << S.completeOrthogonalDecomposition().pseudoInverse().norm() *
    //           S.norm()
    //      << endl;
    // K = _e->E * H.transpose() * S.inverse();
    cout << "S determinant: " << S.sum() << endl;
    cout << " H sum: " << H.array().abs().sum() << endl;
    // K = _e->E * (H.array().colwise() / S.array()).transpose().matrix();
    /*K = _e->E * H.transpose() * S.diagonal().asDiagonal().inverse();
    // K = MatrixXd::Zero(6, m);
    // for (i = 0; i < m; i++)
    //  K += _e->E * H.transpose() / S(i);
    //_e->E = (_e->m->I - K * H) * _e->E;
    //_e->E = _e->E + H*K*J_e->m->_Q_a;
    L = (_e->m->I - K * H); // *_e->E;

    cout << "E Eigenvalues before: " << _e->E.eigenvalues() << endl;
    cout << "KJWJtKt Eigenvalues: "
         << (K * J * _e->m->_W_a * J.transpose() * K.transpose())
              .eigenvalues()
              .transpose()
         << endl;
    cout << "LELt Eigenvalues: "
         << (L * _e->E * L.transpose()).eigenvalues().transpose() << endl;
    _e->E = L * _e->E * L.transpose() +
            K * J * _e->m->_W_a * J.transpose() * K.transpose();*/
    // this kind of "works"
    for (int i = 0; i < H.rows(); i++) {
      K = _e->E * H(i, all).transpose() / S(i, i);
      L = _e->m->I - K * H(i, all);
      _e->E = L * _e->E * L.transpose() + K * J(i, all) * _e->m->_W_a *
                                              J(i, all).transpose() *
                                              K.transpose();
      _e->E = (_e->E + _e->E.transpose()) / 2;
    }
    cout << "E Eigenvalues after: " << _e->E.eigenvalues().transpose() << endl;
    cout << "E: " << endl << _e->E.diagonal().transpose() << endl;
    // cout << "E-VAL: " << _e->E.diagonal().array().abs().sum() << endl;
    if (!isPsd(_e->E)) {
      throw std::runtime_error("Non positive semi-definite matrix!");
    }
    // system propagation
    _e->E = _e->E + _e->m->_Q_a;
    return 0;
  }

  // private:
  // Matrix<double, Dynamic, 12> df;
  // attributes
  int j, i;
  Entity *_e;
  MatrixXd _mk, _m0, _df, S;
  MatrixXd H, J, K, L;
  VectorXd _pose, _pk, _p0, _r, _dr;
  // Matrix<double, Dynamic, 1> S;
  ofstream odbg;
};

#endif // EKF_HPP
