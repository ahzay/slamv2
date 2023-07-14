//
// Created by user on 7/13/23.
//

#include "linemodel.h"
#include <Eigen/Eigenvalues>
#include "entity.h"

LineModel::LineModel(const string &file) : Model(file) {}


Matrix<double, 1, -1> LineModel::dfs(const Entity &e, const Data &d) const {
    const auto measurement = d.get_rotated_measurement();
    ad::var
            r(e._p(0)), a(e._p(1)), xp(d._pose(0)), yp(d._pose(1)),
            dd(measurement(0)), an(measurement(1)), mud(0.),
            muan(0.);
    ad::var f = ft(r, a, xp, yp, dd, an, mud, muan);

    Vector<double, 8> ans;
    auto [dr, da, dxp, dyp, ddd, dan, dmud, dmuan] =
            ad::derivatives(f, ad::wrt(r, a, xp, yp, dd, an, mud, muan));
    ans(0) = dr;
    ans(1) = da;
    ans(2) = dxp;
    ans(3) = dyp;
    ans(4) = ddd;
    ans(5) = dan;
    ans(6) = dmud;
    ans(7) = dmuan;
    return ans;
}

VectorXd LineModel::init(const Aggregate &a) const {
    return Eigen::VectorXd();
}

void LineModel::ls(Entity &e, const Aggregate &a) const {
    Vector<double, 2> p;
    MatrixXd X = a.get_xy_mat();
    MatrixXd b = (X.transpose() * X).inverse() * X.transpose() *
                 Matrix<double, Dynamic, 1>::Ones(X.rows(), 1);
    double r = sgn(b(0)) / sqrt(b(0) * b(0) + b(1) * b(1));
    double aa = atan2(b(1) * r, b(0) * r);
    VectorXd angles = atan2(X.col(1), X.col(0));
    p << r, aa; //, angles.minCoeff(), angles.maxCoeff();
    // cout << p.transpose() << endl;
    e._p = p;
    dop(e,a);
}

struct LineModel::CostFunction {
    CostFunction(double x_p,
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
        residual[0] = ((_x * cos(p[1]) + _y * sin(p[1])) / p[0]) - 1.;
        //    measurement variation constraint
        residual[1] = abs((d[0] - _ap_d) * (d[0] - _ap_d) / _sigma_d);
        // residual[1] *= 100000;
        //   parameter variation constraint
        Eigen::Map<const Eigen::Vector<T, 2>> pv(p);
        Eigen::Vector<T, 2> dp = pv - _ap_p;
        residual[2] = abs((dp.transpose() * _ap_dop.inverse() * dp)(0));
        return true;
    }

private:
    double _ap_d, _an, _x_p, _y_p, _sigma_d;
    VectorXd _ap_p;
    MatrixXd _ap_dop;
};

void LineModel::ap_ls(Entity &e, const Aggregate &a) const {
    cout<<"line ap_ls: "<<endl;
    const auto dan = a.get_rotated_measurement_mat();
    const auto xy = a.get_xy_mat();
    VectorXd d = dan.col(0);
    VectorXd an = dan.col(1);
    VectorXd x = xy.col(0);
    VectorXd y = xy.col(1);
    // TODO: GENERALIZE FOR N-DIM MEASUREMENTS
    ceres::Problem problem;
    // parameter count is the number of model params + number of measurements
    double n_p[2];
    vector<double> n_d(d.size());
    // fill initial values
    for (int i = 0; i < 2; i++)
        n_p[i] = e._p(i);
    for (int i = 0; i < d.size(); i++)
        n_d[i] = d[i];
    // fill problem
    for (unsigned i = 0; i < d.size(); i++) {
        // each parameter block optimizes the 6 parameters and one measurement
        // cout << endl << ap_p << endl << ap_dop << endl;
        ceres::CostFunction *cost_function =
                new ceres::AutoDiffCostFunction<CostFunction, 3, 2, 1>(
                        new CostFunction(a._pose(0), a._pose(1), d(i), an(i),
                                         _dop_sigma, e._p, e._E));
        problem.AddResidualBlock(cost_function, nullptr, n_p, &n_d[i]);
    }
    // bounds
    /*problem.SetParameterLowerBound(n_p, 3, 0.25); // a
    prblem.SetParameterLowerBound(n_p, 4, 0.25); // b
    problem.SetParameterLowerBound(n_p, 5, 0.3);  // eps <- very important
    problem.SetParameterUpperBound(n_p, 3, 30);   // a
    problem.SetParameterUpperBound(n_p, 4, 30);   // b
    problem.SetParameterUpperBound(n_p, 5, 1.7);*/  // eps
    for (int j = 0; j < _parameter_count; j++) {
        problem.SetParameterLowerBound(n_p, j, _parameter_mins[j]);
        problem.SetParameterUpperBound(n_p, j, _parameter_maxs[j]);
    }
    // solve
    ceres::Solver::Options options;
    options.num_threads = 16;
    options.minimizer_progress_to_stdout = false;
    // options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY; //<- slower
    // options.preconditioner_type = ceres::CLUSTER_JACOBI;       //<- also slower
    // options.initial_trust_region_radius = 1e8;                 // important
    options.max_num_iterations = 100;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // fixes
    VectorXd ans = Map<Vector<double, 2>>(n_p);
    //ans.conservativeResize(6 + d.size());
    Aggregate n_a(a);
    for (int i = 0; i < d.size(); i++)
        n_a._data_vector[i]._measurement(0) = n_d[i];
    anglefix:
    if (ans(1) > 2 * M_PI) {
        ans(1) -= 2 * M_PI;
        goto anglefix;
    }
    if (ans(1) < 0) {
        ans(1) += 2 * M_PI;
        goto anglefix;
    }
    cout << "AP_LS YIELDS: " << ans(seq(0, 1)).transpose() << endl;
    e._p = ans(seq(0, 1));
    ap_dop(e, n_a);
}

bool LineModel::safety(Entity &e) const {
    return false;
}

void LineModel::init_post(Entity &e, const Aggregate &a) const {
    auto xy = a.get_xy_mat();
    ArrayXd x = xy.col(0);
    ArrayXd y = xy.col(1);
    // we first project the points(Sophie, eq 5.2)
    ArrayXd x_p =
            x * pow(sin(e._p(1)), 2) - y * sin(e._p(1))
                                       * cos(e._p(1)) + e._p(0) * cos(e._p(1));
    ArrayXd y_p =
            y * pow(cos(e._p(1)), 2) - x * sin(e._p(1)) * cos(e._p(1)) + e._p(0) * sin(e._p(1));
    // we find the angles of all projected points
    ArrayXd angles = atan2(y_p, x_p);
    e._t(0) = angles.colwise().minCoeff()(0);
    e._t(1) = angles.colwise().maxCoeff()(0);
    cout << "min/max is " << e._t.transpose();
}

bool LineModel::associate(const Entity &e, const Data &d) const {
    auto xy = d.get_xy();
    double x = xy(0);
    double y = xy(1);
    // project data
    double x_p =
            x * pow(sin(e._p(1)), 2) - y * sin(e._p(1)) * cos(e._p(1)) + e._p(0) * cos(e._p(1));
    double y_p =
            y * pow(cos(e._p(1)), 2) - x * sin(e._p(1)) * cos(e._p(1)) + e._p(0) * sin(e._p(1));
    double angle = atan2(y_p, x_p);
    double dist_p = sqrt(pow(x_p - x, 2) + pow(y_p - y, 2));
    // double r = sqrt(x_p * x_p + y_p * y_p); // for dist check
    //  check that angle is close enough to min/max or between them
    double angle_tolerance = 0.1;
    double dist_tolerance = 0.1;
    if ((angle > e._t(0) - angle_tolerance && angle < e._t(1) + angle_tolerance) &&
        (dist_p < dist_tolerance))
        return true;
    else
        return false;
}

void LineModel::augment_post(Entity &e, const Data &d) const {
    auto xy = d.get_xy();
    double x = xy(0);
    double y = xy(1);
    // project data
    double x_p =
            x * pow(sin(e._p(1)), 2) - y * sin(e._p(1)) * cos(e._p(1)) + e._p(0) * cos(e._p(1));
    double y_p =
            y * pow(cos(e._p(1)), 2) - x * sin(e._p(1)) * cos(e._p(1)) + e._p(0) * sin(e._p(1));
    double angle = atan2(y_p, x_p);
    if (angle < e._t(0))
        e._t(0) = angle;
    if (angle > e._t(1))
        e._t(1) = angle;
}

double LineModel::fss(const VectorXd &p, const VectorXd &pose, const VectorXd &rotated_measurement) const {
    return ft<double>(p(0), p(1), pose(0), pose(1), rotated_measurement(0),
                      rotated_measurement(1), 0, 0);
}

template<typename T>
T LineModel::ft(T r, T a, T xp, T yp, T d, T an, T mud, T muan) const {
    T x = xp + (d + mud) * cos(an + muan);
    T y = yp + (d + mud) * sin(an + muan);
    /*T testing_angle = atan2(y, x);
    if (testing_angle >= min(th_min, th_max) &&
        testing_angle <= max(th_min, th_max)) {
      return (x * cos(a) + y * sin(a)) / r - 1.;
    } else
      return 0.;*/
    return ((x * cos(a) + y * sin(a)) / r) - 1.;
}

