#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "preproc.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <autodiff/reverse/var.hpp>
#include <ceres/autodiff_cost_function.h>
#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>
// dogshit
/*bool
areVectorsIndependent(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2)
{
  if (v1.isZero(1e-26) || v2.isZero(1e-26)) {
    return false; // A zero vector is not independent
  }

  double ratio = v1[0] / v2[0];
  for (int i = 1; i < v1.size(); ++i) {
    if (std::abs(v1[i] / v2[i] - ratio) > 1e-26) {
      return true; // Vectors are independent
    }
  }

  return false; // Vectors are dependent
}

bool
isIndependentFromSet(const Eigen::VectorXd& v,
                     const std::vector<Eigen::VectorXd>& set)
{
  for (const auto& vec : set) {
    if (!areVectorsIndependent(v, vec)) {
      return false;
    }
  }
  return true;
}

int
rowRank(const Eigen::MatrixXd& matrix)
{
  std::vector<Eigen::VectorXd> independentRows;
  for (int i = 0; i < matrix.rows(); ++i) {
    Eigen::VectorXd row = matrix.row(i);
    if (isIndependentFromSet(row, independentRows)) {
      independentRows.push_back(row);
    }
  }
  return independentRows.size();
}*/
Eigen::MatrixXd
gramSchmidt(const Eigen::MatrixXd &A) {
    Eigen::MatrixXd B = A;
    for (int i = 0; i < B.cols(); ++i) {
        for (int j = 0; j < i; ++j) {
            B.col(i) -=
                    B.col(j) * (B.col(j).dot(B.col(i))) / (B.col(j).dot(B.col(j)));
        }
    }
    return B;
}

int
rowRank(const Eigen::MatrixXd &A) {
    Eigen::MatrixXd B = gramSchmidt(A.transpose());
    int rank = 0;
    for (int i = 0; i < B.cols(); ++i) {
        if (!B.col(i).isZero(1e-24)) {
            rank++;
        }
    }
    return rank;
}

// predeclarations because clusterfuck, separate .h/.cpp
//******* model independant parameters *******//
template<class MatrixT>
bool
isPsd(const MatrixT &A) {
    // if (!A.isApprox(A.transpose(), 1e-8f)) {
    //   cout << "is not symmetric!" << endl;
    //   cout << "A-A':" << endl << A - A.transpose() << endl;
    //   return false;
    // }
    const auto ldlt = A.template selfadjointView<Eigen::Upper>().ldlt();
    if (ldlt.info() == Eigen::NumericalIssue || !ldlt.isPositive()) {
        cout << "eigenvalues: " << endl << A.eigenvalues().transpose() << endl;
        return false;
    }
    return true;
}

// continuity
double angle_tolerance = 0.5, dist_tolerance = 1;
unsigned short N = 15;
// association mahalanobis
double mahalanobis_aug = 1;

//********************************************//
class Entity;

class Model;

// one measurement
class Data {
public:
    Data() {}

    Data(VectorXd data, VectorXd pose) {
        d = data;
        p = pose;
    }

    // attributes
    VectorXd d;          // data (measurement)
    VectorXd p;          // position at measurement
    Entity *e = nullptr; // associated entity
};

// aggregate to be estimated
class Aggregate {
public:
    // functions
    MatrixXd get_mat() const { // returns a matrix of measurements (for ls/dop)
        MatrixXd m;
        m.conservativeResize(data_vector.size(), 4); // x,y,d,an
        for (int i = 0; i < data_vector.size(); i++)
            m.row(i) = data_vector[i].d;
        return m; // optimize this fun?
    }

    void push_back(Data data) { data_vector.push_back(data); }

    bool check_continuity(Data data) {
        double dist_variation = abs(data.d(2) - data_vector.back().d(2));
        double angle_variation = abs(atan2(sin(data.d(3) - data_vector.back().d(3)),
                                           cos(data.d(3) - data_vector.back().d(3))));
        // double angle_variation = abs(v(3) - data(last, 3));
        //  ret 1 for error
        //  PRINT FOR DEBUG
        bool ans =
                angle_variation > angle_tolerance || dist_variation > dist_tolerance;
        if (ans)
            cout << "Dist variation: " << dist_variation << endl
                 << "Angle variation: " << angle_variation << endl
                 << "Dist tolerance " << dist_tolerance << endl
                 << "Angle tolerance " << angle_tolerance << endl;
        return ans;
    }

    void flush(Data data) {                 // flush and fill with one data
        data_vector.clear(); // resetting
        data_vector.push_back(data);
        position = data.p;
        // e = nullptr;
        m = nullptr;
    }

    // attributes
    Model *m = nullptr;
    // Entity *e = nullptr;
    //  MatrixXd data; // does not include position
    VectorXd position; // not necessary
    vector<Data> data_vector;
};

//
// static Visualizer v;
// simple functions
// sign
template<typename T>
int
sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double
cov(VectorXd x, VectorXd y) {
    auto xm = x.array() - x.mean();
    auto ym = y.array() - y.mean();
    auto m = xm.array() * ym.array();
    return m.sum() / (m.size() - 1);
}

VectorXd
atan2(VectorXd y, VectorXd x) {
    VectorXd res(y.size());
    for (int i = 0; i < y.size(); i++) {
        res(i) = atan2(y(i), x(i));
    }
    return res;
}

namespace ad = autodiff;

// MODEL 0: ELLIPSES
template<typename T>
T
ft_0(T xc, T yc, T th, T a, T b, T e, T xp, T yp, T d, T an, T mud, T muan) {
    T x = xp + (d + mud) * cos(an + muan);
    T y = yp + (d + mud) * sin(an + muan);
    T f1 = ((x - xc) * cos(th) + (y - yc) * sin(th)) / a;
    T f2 = ((x - xc) * sin(th) - (y - yc) * cos(th)) / b;
    T buf = pow(pow(f1, 2.), 1. / e) + pow(pow(f2, 2.), 1. / e);
    return sqrt(a * b) * (pow(buf, e) - 1.);
}

template<typename T>
T
fs_0(const VectorXd p, const VectorXd pos, const VectorXd data) {
    // return ft<T, T, T, T, T, T, T, T, T, T, T, T>(xc, yc, th, a, b, e, xp,
    // yp, d, an, 0, 0);
    return ft_0<T>(
            p(0), p(1), p(2), p(3), p(4), p(5), pos(0), pos(1), data(2), data(3), 0, 0);
}

/*! \fn VectorXd fsn_0(const Vector<double, 6> &p, const Vector2d &pos, const
   Matrix<double, Dynamic, 2> &m) \
    \brief Implicit ellipse function
    \param p function parameters
    \param pos position of scanner
    \param m matrix of measurements
    \return returns error
*/
VectorXd
fsn_0(const VectorXd p, const VectorXd pos, const MatrixXd data) {
    VectorXd ans(data.rows());
    for (int i = 0; i < data.rows(); i++)
        ans(i) = fs_0<double>(p, pos, data(i, all));
    return ans;
}

/*! \fn VectorXd dfsn_0(VectorXd p, VectorXd pos, MatrixXd m)
    \brief Implicit ellipse derivated function
    \param p function parameters
    \param pos position of scanner
    \param m matrix of measurements
    \return returns implicit derivative
*/
MatrixXd
dfsn_0(const VectorXd p, const VectorXd pos, const MatrixXd data) {
    ad::var xc(p(0)), yc(p(1)), th(p(2)), a(p(3)), b(p(4)), e(p(5)), xp(pos(0)),
            yp(pos(1)), d(0.), an(0.), mud(0.), muan(0.);
    ad::var f = ft_0(xc, yc, th, a, b, e, xp, yp, d, an, mud, muan);

    MatrixXd ans(data.rows(), 12);
    for (unsigned i = 0; i < data.rows(); i++) {
        d.update(data.row(i)(2));
        an.update(data.row(i)(3));
        f.update();
        auto [dxc, dyc, dth, da, db, de, dxp, dyp, dd, dan, dmud, dmuan] =
                ad::derivatives(f,
                                ad::wrt(xc, yc, th, a, b, e, xp, yp, d, an, mud, muan));
        ans(i, 0) = dxc;
        ans(i, 1) = dyc;
        ans(i, 2) = dth;
        ans(i, 3) = da;
        ans(i, 4) = db;
        ans(i, 5) = de;
        ans(i, 6) = dxp;
        ans(i, 7) = dyp;
        ans(i, 8) = dd;
        ans(i, 9) = dan;
        ans(i, 10) = dmud;
        ans(i, 11) = dmuan;
        // cout << "data row: " << data.row(i) << endl;
        // cout << "params: " << p.transpose() << endl;
        // cout << "position: " << pos.transpose() << endl;
        // cout << "df row: " << ans(i, all) << endl;
    }
    return ans;
}

MatrixXd
dop_0(const VectorXd p,
      const VectorXd pos,
      const MatrixXd data,
      const double sigma_d) {
    Matrix<double, Dynamic, 12> df = dfsn_0(p, pos, data);
    MatrixXd Jmes = df(all, {10});
    // VectorXd Emes = Jmes.array() * sigma_d * sigma_d * Jmes.array();
    MatrixXd H = df(all, {0, 1, 2, 3, 4, 5});
    // cout << "HtH: " << endl << H.transpose() * H << endl;
    // cout << "JtJ: " << endl << Jmes.transpose() * Jmes << endl;

    // MatrixXd En = (H.transpose() *
    //                (Jmes * sigma_d * sigma_d * Jmes.transpose()).inverse() * H)
    //                 .inverse();
    // MatrixXd En = ((H.array().colwise() / Emes.array()).matrix().transpose() *
    //               (H.array().colwise() / Emes.array()).matrix())
    //                .inverse();
    MatrixXd En = MatrixXd::Zero(6, 6);
    for (int i = 0; i < df.rows(); i++)
        En +=
                (H(i, all).transpose() *
                 (Jmes(i, all) * sigma_d * sigma_d * Jmes(i, all).transpose()).inverse() *
                 H(i, all));
    En = En.inverse();
    cout << "En: " << endl << En << endl;
    if (!isPsd(En)) {
        throw std::runtime_error("Non positive semi-definite matrix!");
    }
    return En;
}

VectorXd
init_0(VectorXd x, VectorXd y, VectorXd loc) {
    double *p = new double[11];
    MatrixXd m({{cov(x, x), cov(x, y)},
                {cov(y, x), cov(y, y)}});
    JacobiSVD<MatrixXd> svd(m, ComputeFullU | ComputeFullV);
    auto V = svd.matrixV();
    auto D = svd.singularValues();
    double cond = D(0) / D(D.size() - 1);
    double th0 = atan2(V(1, 0), V(0, 0));
    // rotation of ptcloud
    MatrixXd d(x.size(), 2);
    d.col(0) = x;
    d.col(1) = y;
    MatrixXd rotm(2, 2);
    rotm << cos(-th0), -sin(-th0), sin(-th0), cos(-th0);
    MatrixXd rotd = (rotm * d.transpose()).transpose();
    double x0 = d.col(0).mean();
    double y0 = d.col(1).mean();
    double a0 = (rotd.col(0).maxCoeff() - rotd.col(0).minCoeff()) / 2;
    double b0 = (rotd.col(1).maxCoeff() - rotd.col(1).minCoeff()) / 2;
    double angle = atan2(y0 - loc(1), x0 - loc(0));
    p[0] = x0 + sqrt(a0 + b0) * cos(angle); // pushes too much !!!
    p[1] = y0 + sqrt(a0 + b0) * sin(angle);
    p[6] = x0 + sqrt(a0 * a0 + b0 * b0) * cos(angle); // alt x0 and y0
    p[7] = y0 + sqrt(a0 * a0 + b0 * b0) * sin(angle);
    // using old a0 b0
    a0 = sqrt(2 * D(0)); // a0
    b0 = sqrt(2 * D(1)); // b0
    //
    if (b0 < 0.25)
        b0 = 0.25;
    if (a0 < 0.25)
        a0 = 0.25;
    p[2] = th0;
    p[3] = a0;
    p[4] = b0;
    p[5] = 1.0;            // CHANGED
    p[8] = th0 + M_PI / 2; // alt th0
    p[9] = x0;
    p[10] = y0; // classic x0 and y0
    VectorXd p_i = Map<Vector<double, 11>>(p);
    delete[] p;
    return p_i;
}

struct LossFunction_0 {
    LossFunction_0(double x, double y)
            : _x(x), _y(y) {
    }

    template<typename T>
    bool operator()(const T *const p, T *residual) const {
        T f1 = ((_x - p[0]) * cos(p[2]) + (_y - p[1]) * sin(p[2])) / p[3];
        T f2 = ((_x - p[0]) * sin(p[2]) - (_y - p[1]) * cos(p[2])) / p[4];
        residual[0] = // log(
                pow((1. + p[4]) * (1. + p[3]),
                    2.) // penalty on area
                // exp(abs(p[5] - 1.0)) * // penalise eps, pow 1 = no penalty
                * (pow(pow(pow(f1, 2.0), (1.0 / p[5])) + pow(pow(f2, 2.0), (1.0 / p[5])),
                       p[5]) -
                   1.); // +
        // 1.);
        return true;
    }

private:
    double _x;
    double _y;
};

VectorXd
ls_0(const VectorXd loc, const MatrixXd data) {
    VectorXd xdata = data(all, 0);
    VectorXd ydata = data(all, 1);
    auto p0 = init_0(xdata, ydata, loc);
    double p0arr[6] = {p0[0], p0[1], p0[2], p0[3], p0[4], p0[5]};
    // dist variation
    // double p1[6] = {p0[6], p0[7], p0[2], p0[3], p0[4], p0[5]};
    // theta variation
    double p2[6] = {p0[0], p0[1], p0[8], p0[3], p0[4], p0[5]};
    // double p3[6] = {p1[0], p1[1], p1[8], p1[3], p1[4], p0[5]};
    //
    const int n = 2;
    // double *pa[n] = {p8, p2};
    // double *pa[n] = {p0arr, p1, p2, p3};
    double *pa[n] = {p0arr, p2};
    //   double *pa[] = {p0.data(), p1, p2, p3};
    // double* pa[] = { p2, p4 };
    //  ceres::Problem *problem = new ceres::Problem[n]; // destructor segfaults
    ceres::Problem problem[n];
    for (unsigned i = 0; i < xdata.size(); i++)
        for (unsigned k = 0; k < n; k++) {
            ceres::CostFunction *cost_function =
                    new ceres::AutoDiffCostFunction<LossFunction_0, 1, 6>(
                            new LossFunction_0(xdata(i), ydata(i)));
            problem[k].AddResidualBlock(cost_function, nullptr, pa[k]);
        }
    // bounds
    for (unsigned i = 0; i < n; i++) {
        problem[i].SetParameterLowerBound(pa[i], 3, 0.25); // a
        problem[i].SetParameterLowerBound(pa[i], 4, 0.25); // b
        problem[i].SetParameterLowerBound(pa[i], 5, 0.3);  // eps <- very important
        problem[i].SetParameterUpperBound(pa[i], 3, 30);   // a
        problem[i].SetParameterUpperBound(pa[i], 4, 30);   // b
        problem[i].SetParameterUpperBound(pa[i], 5, 1.7);  // eps
        // avoids auto-occlusion
        if (loc[0] - p0[9] > 0) // robot on right of shape
            problem[i].SetParameterUpperBound(pa[i], 0, p0[9]); // x
        else
            problem[i].SetParameterLowerBound(pa[i], 0, p0[9]);
        if (loc[1] - p0[10] > 0)                               // robot above shape
            problem[i].SetParameterUpperBound(pa[i], 1, p0[10]); // y
        else
            problem[i].SetParameterLowerBound(pa[i], 1, p0[10]);
    }
    ceres::Solver::Options options;
    options.num_threads = 24;
    options.minimizer_progress_to_stdout = false;
    options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY; //<- slower
    options.preconditioner_type = ceres::CLUSTER_JACOBI;       //<- also slower
    options.initial_trust_region_radius = 1e8;                 // important
    options.max_num_iterations = 100;
    ceres::Solver::Summary summary[n];
    for (unsigned i = 0; i < n; i++)
        ceres::Solve(options, &problem[i], &summary[i]);
    unsigned opt = 0;
    double min_cost = 10000;

    for (unsigned i = 0; i < n; i++) {
        cout << pa[i][0] << ", " << pa[i][1] << ", " << pa[i][2] << ", " << pa[i][3]
             << ", " << pa[i][4] << ", " << pa[i][5] << endl;
        // checking that array is finite (not nan or inf)
        bool arrayIsFinite = 1;
        for (int k = 0; k < 6; k++)
            if (!isfinite(pa[i][k])) {
                arrayIsFinite = 0;
                break;
            }
        if (!arrayIsFinite)
            continue;
        // if not skipped in continue, has a chance at being optimal sol
        if (summary[i].final_cost < min_cost && arrayIsFinite) {
            min_cost = summary[i].final_cost;
            opt = i;
        }
    }
    Vector<double, 6> ans = Map<Vector<double, 6>>(pa[opt]);
// fix angle
    anglefix:
    if (ans(2) > 2 * M_PI) {
        ans(2) -= 2 * M_PI;
        goto anglefix;
    }
    if (ans(2) < 0) {
        ans(2) += 2 * M_PI;
        goto anglefix;
    }
    cout << "LS YIELDS: " << ans.transpose() << endl;
    return ans;
}

struct Cost_0 {
    Cost_0(double x_p,
           double y_p,
           double ap_d,
           double an,
           double sigma_d,
           VectorXd ap_p,
           MatrixXd ap_dop)
            : _ap_d(ap_d), _an(an), _x_p(x_p), _y_p(y_p), _sigma_d(sigma_d), _ap_p(ap_p), _ap_dop(ap_dop) {
    }

    template<typename T>
    bool operator()(const T *const p, const T *const d, T *residual) const {
        // nonlin eq constraint
        T _x = _x_p + d[0] * cos(_an);
        T _y = _y_p + d[0] * sin(_an);
        T f1 = ((_x - p[0]) * cos(p[2]) + (_y - p[1]) * sin(p[2])) / p[3];
        T f2 = ((_x - p[0]) * sin(p[2]) - (_y - p[1]) * cos(p[2])) / p[4];
        residual[0] = // log(
                // penalty on area
                abs(pow((1. + p[4]) * (1. + p[3]), 0.5)) *
                //* pow(2, abs(p[5] - 1.0)) * // penalise eps, pow 1 = no
                //   penalty
                abs(
                        (pow(pow(pow(f1, 2.0), (1.0 / p[5])) + pow(pow(f2, 2.0), (1.0 / p[5])),
                             p[5]) -
                         1.)); // +
        // 1.);
        residual[0] *= 1.;
        //    measurement variation constraint
        residual[1] = abs((d[0] - _ap_d) * (d[0] - _ap_d) / _sigma_d);
        // residual[1] *= 100000;
        //   parameter variation constraint
        Eigen::Map<const Eigen::Vector<T, 6>> pv(p);
        Eigen::Vector<T, 6> dp = pv - _ap_p;
        residual[2] = abs((dp.transpose() * _ap_dop.inverse() * dp)(0));
        // residual[2] *= 100000;
        //   area
        residual[3] = 0. * abs(pow((1. + p[4]) * (1. + p[3]), 0.1));
        // residual[3] = 0.;
        return true;
    }

private:
    double _ap_d, _an, _x_p, _y_p, _sigma_d;
    VectorXd _ap_p;
    MatrixXd _ap_dop;
};

VectorXd
ap_ls_0(const VectorXd ap_p,   // apriori params
        const MatrixXd ap_dop, // apriori dop
        const VectorXd loc,
        const MatrixXd xydan,
        const double sigma_d) // measurements
{

    VectorXd d = xydan(all, 2);
    VectorXd an = xydan(all, 3);
    VectorXd x = loc[0] + d.array() * an.array().cos();
    VectorXd y = loc[1] + d.array() * an.array().sin();
    // TODO: GENERALIZE FOR N-DIM MEASUREMENTS
    ceres::Problem problem;
    // parameter count is the number of model params + number of measurements
    double n_p[6];
    vector<double> n_d(d.size());
    // fill initial values
    for (int i = 0; i < 6; i++)
        n_p[i] = ap_p[i];
    for (int i = 0; i < d.size(); i++)
        n_d[i] = d[i];
    // fill problem
    for (unsigned i = 0; i < d.size(); i++) {
        // each parameter block optimizes the 6 parameters and one measurement
        // cout << endl << ap_p << endl << ap_dop << endl;
        ceres::CostFunction *cost_function =
                new ceres::AutoDiffCostFunction<Cost_0, 4, 6, 1>(
                        new Cost_0(loc(0), loc(1), d(i), an(i), sigma_d, ap_p, ap_dop));
        problem.AddResidualBlock(cost_function, nullptr, n_p, &n_d[i]);
    }
    // bounds
    problem.SetParameterLowerBound(n_p, 3, 0.25); // a
    problem.SetParameterLowerBound(n_p, 4, 0.25); // b
    problem.SetParameterLowerBound(n_p, 5, 0.3);  // eps <- very important
    problem.SetParameterUpperBound(n_p, 3, 30);   // a
    problem.SetParameterUpperBound(n_p, 4, 30);   // b
    problem.SetParameterUpperBound(n_p, 5, 1.7);  // eps
    // avoids auto-occlusion, only if enough points
    // if enough points we assume that we see enough of the ellipse
    // for an accurate mean
    // TODO: change to N
    if (x.size() > 35) {
        if (loc[0] - x.mean() > 0) // robot on right of shape
            problem.SetParameterUpperBound(n_p, 0, x.mean()); // x
        else
            problem.SetParameterLowerBound(n_p, 0, x.mean());
        if (loc[1] - y.mean() > 0)                          // robot above shape
            problem.SetParameterUpperBound(n_p, 1, y.mean()); // y
        else
            problem.SetParameterLowerBound(n_p, 1, y.mean());
    }
    // solve
    ceres::Solver::Options options;
    options.num_threads = 24;
    options.minimizer_progress_to_stdout = false;
    // options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY; //<- slower
    // options.preconditioner_type = ceres::CLUSTER_JACOBI;       //<- also slower
    // options.initial_trust_region_radius = 1e8;                 // important
    options.max_num_iterations = 100;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // fixes
    VectorXd ans = Map<Vector<double, 6>>(n_p);
    ans.conservativeResize(6 + d.size());
    for (int i = 0; i < d.size(); i++)
        ans[i + 6] = n_d[i];
    anglefix:
    if (ans(2) > 2 * M_PI) {
        ans(2) -= 2 * M_PI;
        goto anglefix;
    }
    if (ans(2) < 0) {
        ans(2) += 2 * M_PI;
        goto anglefix;
    }
    cout << "AP_LS YIELDS: " << ans(seq(0, 5)).transpose() << endl;
    return ans;
    // new DOP has to be calculated externally
}

bool
safety_0(VectorXd &p) {
    if (p(5) < 0.3) {
        p(5) = 0.3;
        return true;
    }
    if (p(5) > 1.7) {
        p(5) = 1.7;
        return true;
    }
    if (p(3) < 0.25) {
        p(3) = 0.25;
        return true;
    }
    if (p(4) < 0.25) {
        p(4) = 0.25;
        return true;
    }
    if (p(2) > 2 * M_PI) {
        p(2) -= 2 * M_PI;
        return false;
    }
    if (p(2) < 0) {
        p(2) += 2 * M_PI;
        return false;
    }
    return false;
}

VectorXd
dst_0(const VectorXd p, const VectorXd pos, const MatrixXd data) {
    VectorXd ans(data.rows());
    for (int i = 0; i < data.rows(); i++) {
        VectorXd row = data(i, all);
        Vector2d pt = row({0, 1}); // point to find distance with
        pt = pt - p({0, 1});
        // undo theta rotation
        Matrix<double, 2, 2> rotm(
                {{cos(-p[2]), -sin(-p[2])},
                 {sin(-p[2]), cos(-p[2])}});
        pt = rotm * pt;
        // (13) in Rosin and West
        // assign parameters for clarity
        double a = p(3), b = p(4), e = p(5), xp = pt(0), yp = pt(1);
        /*double ddeno1 = 1. / pow(a, 2. / e);
      double ddeno2 = pow(abs(yp / (xp * b)), 2. / e);
      double deno = ddeno1 + ddeno2;
      double xs = pow(abs(1. / deno), e / 2.);
      double ys = xs * yp / xp;*/

        // my own
        // find angle of data to center of ellipse
        double an = atan2(yp, xp);
        double xs = pow(abs(cos(an)), e) * a * sgn(cos(an));
        double ys = pow(abs(sin(an)), e) * b * sgn(sin(an));

        Vector2d pt_s = {xs, ys}; // point on the ellipse curve
        double dans = (pt_s - pt).norm();
        ans(i) = dans;
    }
    return ans;
}

void
init_post_0(const VectorXd p, VectorXd &t, const Aggregate a) {
}

bool
associate_0(const VectorXd p, VectorXd &t, const Data d) {
    return 1;
}

void
augment_post_0(const VectorXd p, VectorXd &t, const Data d) {
}

// MODEL 1: LINES (not passing through origin)
template<typename T>
// (p.71) p_seg][r,a]
T
ft_1(T r, T a, T xp, T yp, T d, T an, T mud, T muan) {
    T x = xp + (d + mud) * cos(an + muan);
    T y = yp + (d + mud) * sin(an + muan);
    T testing_angle = atan2(y, x);
    /*if (testing_angle >= min(th_min, th_max) &&
        testing_angle <= max(th_min, th_max)) {
      return (x * cos(a) + y * sin(a)) / r - 1.;
    } else
      return 0.;*/
    return ((x * cos(a) + y * sin(a)) / r) - 1.;
}

template<typename T>
T
fs_1(const VectorXd p, const VectorXd pos, const VectorXd data) {
    // return ft<T, T, T, T, T, T, T, T, T, T, T, T>(xc, yc, th, a, b, e, xp,
    // yp, d,
    //                                              an, 0, 0);
    return ft_1<T>(p(0), p(1), pos(0), pos(1), data(2), data(3), 0, 0);
}

/*! \fn VectorXd fsn_0(const Vector<double, 6> &p, const Vector2d &pos, const
   Matrix<double, Dynamic, 2> &m) \
    \brief Implicit ellipse function
    \param p function parameters
    \param pos position of scanner
    \param m matrix of measurements
    \return returns error
*/
VectorXd
fsn_1(const VectorXd p, const VectorXd pos, const MatrixXd data) {
    VectorXd ans(data.rows());
    for (int i = 0; i < data.rows(); i++)
        ans(i) = fs_1<double>(p, pos, data(i, all));
    return ans;
}

/*! \fn VectorXd dfsn_0(VectorXd p, VectorXd pos, MatrixXd m)
    \brief Implicit ellipse derivated function
    \param p function parameters
    \param pos position of scanner
    \param m matrix of measurements
    \return returns implicit derivative
*/
MatrixXd
dfsn_1(const VectorXd p, const VectorXd pos, const MatrixXd data) {
    ad::var r(p(0)), a(p(1)), xp(pos(0)), yp(pos(1)), d(0.), an(0.), mud(0.),
            muan(0.);
    ad::var f = ft_1(r, a, xp, yp, d, an, mud, muan);

    MatrixXd ans(data.rows(), 8);
    for (unsigned i = 0; i < data.rows(); i++) {
        d.update(data.row(i)(2));
        an.update(data.row(i)(3));
        f.update();
        auto [dr, da, dxp, dyp, dd, dan, dmud, dmuan] =
                ad::derivatives(f, ad::wrt(r, a, xp, yp, d, an, mud, muan));
        ans(i, 0) = dr;
        ans(i, 1) = da;
        ans(i, 2) = dxp;
        ans(i, 3) = dyp;
        ans(i, 4) = dd;
        ans(i, 5) = dan;
        ans(i, 6) = dmud;
        ans(i, 7) = dmuan;
    }
    return ans;
}

MatrixXd
dop_1(const VectorXd p,
      const VectorXd pos,
      const MatrixXd data,
      const double sigma_d) {
    Matrix<double, Dynamic, 8> df = dfsn_1(p, pos, data);
    MatrixXd Jmes = df(all, {6, 7});
    MatrixXd Emes = Jmes * sigma_d * sigma_d * Jmes.transpose();
    MatrixXd H = df(all, {0, 1});
    MatrixXd En =
            (H.transpose() * Emes.completeOrthogonalDecomposition().pseudoInverse() * H)
                    .completeOrthogonalDecomposition()
                    .pseudoInverse();
    return En;
}

VectorXd
ls_1(const VectorXd loc, const MatrixXd data) {
    Vector<double, 2> p;
    MatrixXd X = data(all, {0, 1});
    MatrixXd b = (X.transpose() * X).inverse() * X.transpose() *
                 Matrix<double, Dynamic, 1>::Ones(X.rows(), 1);
    double r = sgn(b(0)) / sqrt(b(0) * b(0) + b(1) * b(1));
    double a = atan2(b(1) * r, b(0) * r);
    VectorXd angles = atan2(X.col(1), X.col(0));
    p << r, a; //, angles.minCoeff(), angles.maxCoeff();
    // cout << p.transpose() << endl;
    return p;
}

bool
safety_1(VectorXd &p) {
    return 0;
}

void
init_post_1(const VectorXd p, VectorXd &t, const Aggregate a) {
    ArrayXd x = a.get_mat().col(0);
    ArrayXd y = a.get_mat().col(1);
    // we first project the points(Sophie, eq 5.2)
    ArrayXd x_p =
            x * pow(sin(p(1)), 2) - y * sin(p(1)) * cos(p(1)) + p(0) * cos(p(1));
    ArrayXd y_p =
            y * pow(cos(p(1)), 2) - x * sin(p(1)) * cos(p(1)) + p(0) * sin(p(1));
    // we find the angles of all projected points
    ArrayXd angles = atan2(y_p, x_p);
    t(0) = angles.colwise().minCoeff()(0);
    t(1) = angles.colwise().maxCoeff()(0);
    cout << "min/max is " << t.transpose();
}

bool
associate_1(const VectorXd p, VectorXd &t, const Data d) {
    double x = d.d(0);
    double y = d.d(1);
    // project data
    double x_p =
            x * pow(sin(p(1)), 2) - y * sin(p(1)) * cos(p(1)) + p(0) * cos(p(1));
    double y_p =
            y * pow(cos(p(1)), 2) - x * sin(p(1)) * cos(p(1)) + p(0) * sin(p(1));
    double angle = atan2(y_p, x_p);
    double dist_p = sqrt(pow(x_p - x, 2) + pow(y_p - y, 2));
    // double r = sqrt(x_p * x_p + y_p * y_p); // for dist check
    //  check that angle is close enough to min/max or between them
    if ((angle > t(0) - angle_tolerance && angle < t(1) + angle_tolerance) &&
        (dist_p < dist_tolerance))
        return 1;
    else
        return 0;
}

void
augment_post_1(const VectorXd p, VectorXd &t, const Data d) {
    double x = d.d(0);
    double y = d.d(1);
    // project data
    double x_p =
            x * pow(sin(p(1)), 2) - y * sin(p(1)) * cos(p(1)) + p(0) * cos(p(1));
    double y_p =
            y * pow(cos(p(1)), 2) - x * sin(p(1)) * cos(p(1)) + p(0) * sin(p(1));
    double angle = atan2(y_p, x_p);
    if (angle < t(0))
        t(0) = angle;
    if (angle > t(1))
        t(1) = angle;
}

/**********************************************/

VectorXd
ls_2(const VectorXd loc, const MatrixXd data) {
    VectorXd xdata = data(all, 0);
    VectorXd ydata = data(all, 1);
    auto p0 = init_0(xdata, ydata, loc);
    double p0arr[6] = {p0[0], p0[1], p0[2], p0[3], p0[4], p0[5]};
    // dist variation
    double p1[6] = {p0[6], p0[7], p0[2], p0[3], p0[4], p0[5]};
    // theta variation
    double p2[6] = {p0[0], p0[1], p0[8], p0[3], p0[4], p0[5]};
    double p3[6] = {p1[0], p1[1], p1[8], p1[3], p1[4], p0[5]};
    //
    const int n = 4;
    // double *pa[n] = {p8, p2};
    double *pa[n] = {p0arr, p1, p2, p3};
    //   double *pa[] = {p0.data(), p1, p2, p3};
    // double* pa[] = { p2, p4 };
    //  ceres::Problem *problem = new ceres::Problem[n]; // destructor segfaults
    ceres::Problem problem[n];
    for (unsigned i = 0; i < xdata.size(); i++)
        for (unsigned k = 0; k < n; k++) {
            ceres::CostFunction *cost_function =
                    new ceres::AutoDiffCostFunction<LossFunction_0, 1, 6>(
                            new LossFunction_0(xdata(i), ydata(i)));
            problem[k].AddResidualBlock(cost_function, nullptr, pa[k]);
        }
    // bounds
    for (unsigned i = 0; i < n; i++) {
        problem[i].SetParameterLowerBound(pa[i], 3, 0.25); // a
        problem[i].SetParameterLowerBound(pa[i], 4, 0.25); // b
        problem[i].SetParameterLowerBound(pa[i], 5, 0.8);  // eps <- very important
        problem[i].SetParameterUpperBound(pa[i], 3, 30);   // a
        problem[i].SetParameterUpperBound(pa[i], 4, 30);   // b
        problem[i].SetParameterUpperBound(pa[i], 5, 1.2);  // eps
        // avoids auto-occlusion
        // problem[i].SetParameterUpperBound(pa[i], 0, max(p0[9], p0[6]));  // x
        // problem[i].SetParameterLowerBound(pa[i], 0, min(p0[9], p0[6]));  // x
        // problem[i].SetParameterUpperBound(pa[i], 1, max(p0[10], p0[7])); // y
        // problem[i].SetParameterLowerBound(pa[i], 1, min(p0[10], p0[7])); // y
    }
    ceres::Solver::Options options;
    options.num_threads = 16;
    options.minimizer_progress_to_stdout = false;
    options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY; //<- slower
    options.preconditioner_type = ceres::CLUSTER_JACOBI;       //<- also slower
    options.initial_trust_region_radius = 1e8;                 // important
    options.max_num_iterations = 100;
    ceres::Solver::Summary summary[n];
    for (unsigned i = 0; i < n; i++)
        ceres::Solve(options, &problem[i], &summary[i]);
    unsigned opt = 0;
    double min_cost = 10000;

    for (unsigned i = 0; i < n; i++) {
        cout << pa[i][0] << ", " << pa[i][1] << ", " << pa[i][2] << ", " << pa[i][3]
             << ", " << pa[i][4] << ", " << pa[i][5] << endl;
        // checking that array is finite (not nan or inf)
        bool arrayIsFinite = 1;
        for (int k = 0; k < 6; k++)
            if (!isfinite(pa[i][k])) {
                arrayIsFinite = 0;
                break;
            }
        if (!arrayIsFinite)
            continue;
        // if not skipped in continue, has a chance at being optimal sol
        if (summary[i].final_cost < min_cost && arrayIsFinite) {
            min_cost = summary[i].final_cost;
            opt = i;
        }
    }
    Vector<double, 6> ans = Map<Vector<double, 6>>(pa[opt]);
// fix angle
    anglefix:
    if (ans(2) > 2 * M_PI) {
        ans(2) -= 2 * M_PI;
        goto anglefix;
    }
    if (ans(2) < 0) {
        ans(2) += 2 * M_PI;
        goto anglefix;
    }
    cout << "LS YIELDS: " << ans.transpose() << endl;
    return ans;
}

VectorXd
ap_ls_2(const VectorXd ap_p,   // apriori params
        const MatrixXd ap_dop, // apriori dop
        const VectorXd loc,
        const MatrixXd xydan,
        const double sigma_d) // measurements
{
    VectorXd d = xydan(all, 2);
    VectorXd an = xydan(all, 3);
    // TODO: GENERALIZE FOR N-DIM MEASUREMENTS
    ceres::Problem problem;
    // parameter count is the number of model params + number of measurements
    double n_p[6];
    vector<double> n_d(d.size());
    // fill initial values
    for (int i = 0; i < 6; i++)
        n_p[i] = ap_p[i];
    for (int i = 0; i < d.size(); i++)
        n_d[i] = d[i];
    // fill problem
    for (unsigned i = 0; i < d.size(); i++) {
        // each parameter block optimizes the 6 parameters and one measurement
        // cout << endl << ap_p << endl << ap_dop << endl;
        ceres::CostFunction *cost_function =
                new ceres::AutoDiffCostFunction<Cost_0, 4, 6, 1>(
                        new Cost_0(loc(0), loc(1), d(i), an(i), sigma_d, ap_p, ap_dop));
        problem.AddResidualBlock(cost_function, nullptr, n_p, &n_d[i]);
    }
    // bounds
    problem.SetParameterLowerBound(n_p, 3, 0.25); // a
    problem.SetParameterLowerBound(n_p, 4, 0.25); // b
    problem.SetParameterLowerBound(n_p, 5, 0.8);  // eps <- very important
    problem.SetParameterUpperBound(n_p, 3, 30);   // a
    problem.SetParameterUpperBound(n_p, 4, 30);   // b
    problem.SetParameterUpperBound(n_p, 5, 1.2);  // eps
    // solve
    ceres::Solver::Options options;
    options.num_threads = 16;
    options.minimizer_progress_to_stdout = false;
    // options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY; //<- slower
    // options.preconditioner_type = ceres::CLUSTER_JACOBI;       //<- also slower
    // options.initial_trust_region_radius = 1e8;                 // important
    options.max_num_iterations = 10;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // fixes
    VectorXd ans = Map<Vector<double, 6>>(n_p);
    ans.conservativeResize(6 + d.size());
    for (int i = 0; i < d.size(); i++)
        ans[i + 6] = n_d[i];
    anglefix:
    if (ans(2) > 2 * M_PI) {
        ans(2) -= 2 * M_PI;
        goto anglefix;
    }
    if (ans(2) < 0) {
        ans(2) += 2 * M_PI;
        goto anglefix;
    }
    cout << "AP_LS YIELDS: " << ans(seq(0, 5)).transpose() << endl;
    return ans;
    // new DOP has to be calculated externally
}

bool
safety_2(VectorXd &p) {
    if (p(5) < 0.8) {
        p(5) = 0.8;
        return 1;
    }
    if (p(5) > 1.2) {
        p(5) = 1.2;
        return 1;
    }
    if (p(3) < 0.2) {
        p(3) = 0.2;
        return 1;
    }
    if (p(4) < 0.2) {
        p(4) = 0.2;
        return 1;
    }
    if (p(2) > 2 * M_PI) {
        p(2) -= 2 * M_PI;
        return 0;
    }
    if (p(2) < 0) {
        p(2) += 2 * M_PI;
        return 0;
    }
    return 0;
}

/**********************************************/

VectorXd
ls_3(const VectorXd loc, const MatrixXd data) {
    VectorXd xdata = data(all, 0);
    VectorXd ydata = data(all, 1);
    auto p0 = init_0(xdata, ydata, loc);
    double p0arr[6] = {p0[0], p0[1], p0[2], p0[3], p0[4], p0[5]};
    // dist variation
    double p1[6] = {p0[6], p0[7], p0[2], p0[3], p0[4], p0[5]};
    // theta variation
    double p2[6] = {p0[0], p0[1], p0[8], p0[3], p0[4], p0[5]};
    double p3[6] = {p1[0], p1[1], p1[8], p1[3], p1[4], p0[5]};
    //
    const int n = 4;
    // double *pa[n] = {p8, p2};
    double *pa[n] = {p0arr, p1, p2, p3};
    //   double *pa[] = {p0.data(), p1, p2, p3};
    // double* pa[] = { p2, p4 };
    //  ceres::Problem *problem = new ceres::Problem[n]; // destructor segfaults
    ceres::Problem problem[n];
    for (unsigned i = 0; i < xdata.size(); i++)
        for (unsigned k = 0; k < n; k++) {
            ceres::CostFunction *cost_function =
                    new ceres::AutoDiffCostFunction<LossFunction_0, 1, 6>(
                            new LossFunction_0(xdata(i), ydata(i)));
            problem[k].AddResidualBlock(cost_function, nullptr, pa[k]);
        }
    // bounds
    for (unsigned i = 0; i < n; i++) {
        problem[i].SetParameterLowerBound(pa[i], 3, 0.25); // a
        problem[i].SetParameterLowerBound(pa[i], 4, 0.25); // b
        problem[i].SetParameterLowerBound(pa[i], 5, 1.3);  // eps <- very important
        problem[i].SetParameterUpperBound(pa[i], 3, 30);   // a
        problem[i].SetParameterUpperBound(pa[i], 4, 30);   // b
        problem[i].SetParameterUpperBound(pa[i], 5, 1.7);  // eps
        // avoids auto-occlusion
        // problem[i].SetParameterUpperBound(pa[i], 0, max(p0[9], p0[6]));  // x
        // problem[i].SetParameterLowerBound(pa[i], 0, min(p0[9], p0[6]));  // x
        // problem[i].SetParameterUpperBound(pa[i], 1, max(p0[10], p0[7])); // y
        // problem[i].SetParameterLowerBound(pa[i], 1, min(p0[10], p0[7])); // y
    }
    ceres::Solver::Options options;
    options.num_threads = 16;
    options.minimizer_progress_to_stdout = false;
    options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY; //<- slower
    options.preconditioner_type = ceres::CLUSTER_JACOBI;       //<- also slower
    options.initial_trust_region_radius = 1e8;                 // important
    options.max_num_iterations = 100;
    ceres::Solver::Summary summary[n];
    for (unsigned i = 0; i < n; i++)
        ceres::Solve(options, &problem[i], &summary[i]);
    unsigned opt = 0;
    double min_cost = 10000;

    for (unsigned i = 0; i < n; i++) {
        cout << pa[i][0] << ", " << pa[i][1] << ", " << pa[i][2] << ", " << pa[i][3]
             << ", " << pa[i][4] << ", " << pa[i][5] << endl;
        // checking that array is finite (not nan or inf)
        bool arrayIsFinite = 1;
        for (int k = 0; k < 6; k++)
            if (!isfinite(pa[i][k])) {
                arrayIsFinite = 0;
                break;
            }
        if (!arrayIsFinite)
            continue;
        // if not skipped in continue, has a chance at being optimal sol
        if (summary[i].final_cost < min_cost && arrayIsFinite) {
            min_cost = summary[i].final_cost;
            opt = i;
        }
    }
    Vector<double, 6> ans = Map<Vector<double, 6>>(pa[opt]);
// fix angle
    anglefix:
    if (ans(2) > 2 * M_PI) {
        ans(2) -= 2 * M_PI;
        goto anglefix;
    }
    if (ans(2) < 0) {
        ans(2) += 2 * M_PI;
        goto anglefix;
    }
    cout << "LS YIELDS: " << ans.transpose() << endl;
    return ans;
}

VectorXd
ap_ls_3(const VectorXd ap_p,   // apriori params
        const MatrixXd ap_dop, // apriori dop
        const VectorXd loc,
        const MatrixXd xydan,
        const double sigma_d) // measurements
{
    VectorXd d = xydan(all, 2);
    VectorXd an = xydan(all, 3);
    // TODO: GENERALIZE FOR N-DIM MEASUREMENTS
    ceres::Problem problem;
    // parameter count is the number of model params + number of measurements
    double n_p[6];
    vector<double> n_d(d.size());
    // fill initial values
    for (int i = 0; i < 6; i++)
        n_p[i] = ap_p[i];
    for (int i = 0; i < d.size(); i++)
        n_d[i] = d[i];
    // fill problem
    for (unsigned i = 0; i < d.size(); i++) {
        // each parameter block optimizes the 6 parameters and one measurement
        // cout << endl << ap_p << endl << ap_dop << endl;
        ceres::CostFunction *cost_function =
                new ceres::AutoDiffCostFunction<Cost_0, 4, 6, 1>(
                        new Cost_0(loc(0), loc(1), d(i), an(i), sigma_d, ap_p, ap_dop));
        problem.AddResidualBlock(cost_function, nullptr, n_p, &n_d[i]);
    }
    // bounds
    problem.SetParameterLowerBound(n_p, 3, 0.25); // a
    problem.SetParameterLowerBound(n_p, 4, 0.25); // b
    problem.SetParameterLowerBound(n_p, 5, 1.3);  // eps <- very important
    problem.SetParameterUpperBound(n_p, 3, 30);   // a
    problem.SetParameterUpperBound(n_p, 4, 30);   // b
    problem.SetParameterUpperBound(n_p, 5, 1.7);  // eps
    // solve
    ceres::Solver::Options options;
    options.num_threads = 16;
    options.minimizer_progress_to_stdout = false;
    // options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY; //<- slower
    // options.preconditioner_type = ceres::CLUSTER_JACOBI;       //<- also slower
    // options.initial_trust_region_radius = 1e8;                 // important
    options.max_num_iterations = 10;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // fixes
    VectorXd ans = Map<Vector<double, 6>>(n_p);
    ans.conservativeResize(6 + d.size());
    for (int i = 0; i < d.size(); i++)
        ans[i + 6] = n_d[i];
    anglefix:
    if (ans(2) > 2 * M_PI) {
        ans(2) -= 2 * M_PI;
        goto anglefix;
    }
    if (ans(2) < 0) {
        ans(2) += 2 * M_PI;
        goto anglefix;
    }
    cout << "AP_LS YIELDS: " << ans(seq(0, 5)).transpose() << endl;
    return ans;
    // new DOP has to be calculated externally
}

bool
safety_3(VectorXd &p) {
    if (p(5) < 1.3) {
        p(5) = 1.3;
        return 1;
    }
    if (p(5) > 1.7) {
        p(5) = 1.7;
        return 1;
    }
    if (p(3) < 0.2) {
        p(3) = 0.2;
        return 1;
    }
    if (p(4) < 0.2) {
        p(4) = 0.2;
        return 1;
    }
    if (p(2) > 2 * M_PI) {
        p(2) -= 2 * M_PI;
        return 0;
    }
    if (p(2) < 0) {
        p(2) += 2 * M_PI;
        return 0;
    }
    return 0;
}

#endif // FUNCTIONS_HPP
