//
// Created by user on 7/10/23.
//

#include "ellipsemodel.h"
#include <Eigen/Eigenvalues>
#include "entity.h"


EllipseModel::EllipseModel(const string &file, const CmdLineOptions &options) : Model(file, options) {}

Matrix<double, 1, -1> EllipseModel::dfs(const Entity &e, const Data &d) const {
    const auto measurement = d.get_rotated_measurement();
    ad::var xc(e._p(0)), yc(e._p(1)), th(e._p(2)), aa(e._p(3)), b(e._p(4)), ee(e._p(5)), xp(d._pose(0)),
            yp(d._pose(1)),
            dd(measurement(0)), an(measurement(1)),
            mud(0.), muan(0.);
    ad::var f = ft(xc, yc, th, aa, b, ee, xp, yp, dd, an, mud, muan);
    Vector<double, 12> ans;
    auto [dxc, dyc, dth, da, db, de,
            dxp, dyp, ddd, dan, dmud, dmuan] =
            ad::derivatives(f,
                            ad::wrt(xc, yc, th, aa, b, ee,
                                    xp, yp, dd, an, mud, muan));
    ans(0) = dxc;
    ans(1) = dyc;
    ans(2) = dth;
    ans(3) = da;
    ans(4) = db;
    ans(5) = de;
    ans(6) = dxp;
    ans(7) = dyp;
    ans(8) = ddd;
    ans(9) = dan;
    ans(10) = dmud;
    ans(11) = dmuan;
    return ans;
}


VectorXd EllipseModel::init(const Aggregate &a) const {
    const auto xy = a.get_xy_mat();
    const auto loc = a._pose;
    auto x = xy.col(0), y = xy.col(1);
    auto *p = new double[11];
    MatrixXd m({{vector_cov(x, x), vector_cov(x, y)},
                {vector_cov(y, x), vector_cov(y, y)}});
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
    //p[5] = _parameter_maxs[5];            // CHANGED
    p[5] = (_parameter_mins[5] + _parameter_maxs[5]) / 2;
    p[8] = th0 + M_PI / 2; // alt th0
    p[9] = x0;
    p[10] = y0; // classic x0 and y0
    VectorXd p_i = Map<Vector<double, 11>>(p);
    delete[] p;
    return p_i;
}

struct EllipseModel::LossFunction {
    LossFunction(double x, double y)
            : _x(x), _y(y) {
    }

    template<typename T>
    bool operator()(const T *const p, T *residual) const {
        T f1 = ((_x - p[0]) * cos(p[2]) + (_y - p[1]) * sin(p[2])) / p[3];
        T f2 = ((_x - p[0]) * sin(p[2]) - (_y - p[1]) * cos(p[2])) / p[4];
        residual[0] = 0. *// log(
                      pow((p[4]) * (p[3]), 0.5); // penalty on area (THIS WORKS HORRIBLE)
        residual[1] =
                // exp(abs(p[5] - 1.0)) * // penalise eps, pow 1 = no penalty
                pow((p[4]) * (p[3]), 1.0) *
                (pow(pow(pow(f1, 2.0), (1.0 / p[5])) + pow(pow(f2, 2.0), (1.0 / p[5])),
                     p[5] / 1.) -
                 1.); // +
        //residual[2] = abs(p[3] / p[4] - 1.); // retarded, maybe try with lesser weight?
        // 1.);
        return true;
    }

private:
    double _x;
    double _y;
};

void EllipseModel::ls(Entity &e, const Aggregate &a, const bool &alreadyInitialized) const {
    cout << "starting ellipse ls" << endl;
    const auto xydata = a.get_xy_mat();
    VectorXd xdata = xydata(all, 0);
    VectorXd ydata = xydata(all, 1);
    VectorXd p0 = init(a);
    VectorXd p0o;
    int n;
    const int l = 12; // normal/pushed(p[0,1]) * 0.1/1.0/1.9
    if (alreadyInitialized && 0) {
        p0o = e._p;
        n = l + 1;
    } else {
        p0o = p0;
        n = l;
    }
    // eps var
    double p0arr[6] = {p0[0], p0[1], p0[2], p0[3], p0[4], 1.0};
    double p1arr[6] = {p0[0], p0[1], p0[2], p0[3], p0[4], 0.1};
    double p2arr[6] = {p0[0], p0[1], p0[2], p0[3], p0[4], 1.9};
    // push var
    double p3arr[6] = {p0[6], p0[7], p0[2], p0[3], p0[4], 1.0};
    double p4arr[6] = {p0[6], p0[7], p0[2], p0[3], p0[4], 0.1};
    double p5arr[6] = {p0[6], p0[7], p0[2], p0[3], p0[4], 1.9};
    // rotate var 82% ->
    double p6arr[6] = {p0arr[0], p0arr[1], p0arr[2] + M_PI / 2, p0arr[3], p0arr[4], p0arr[5]};
    double p7arr[6] = {p1arr[0], p1arr[1], p1arr[2] + M_PI / 2, p1arr[3], p1arr[4], p1arr[5]};
    double p8arr[6] = {p2arr[0], p2arr[1], p2arr[2] + M_PI / 2, p2arr[3], p2arr[4], p2arr[5]};
    double p9arr[6] = {p3arr[0], p3arr[1], p3arr[2] + M_PI / 2, p3arr[3], p3arr[4], p3arr[5]};
    double p10arr[6] = {p4arr[0], p4arr[1], p4arr[2] + M_PI / 2, p4arr[3], p4arr[4], p4arr[5]};
    double p11arr[6] = {p5arr[0], p5arr[1], p5arr[2] + M_PI / 2, p5arr[3], p5arr[4], p5arr[5]};
    // alreadyInit
    double p0oarr[6] = {xdata.mean(), ydata.mean(), p0o[2], p0o[3], p0o[4], p0o[5]};
    //double *pa[l+1] = {p0arr,p1arr,p2arr,p3arr,p4arr,p5arr, p0oarr};
    double *pa[l + 1] = {p0arr, p1arr, p2arr, p3arr, p4arr, p5arr, p6arr, p7arr, p8arr, p9arr, p10arr, p11arr, p0oarr};
    // theta variation
    /*double p1[6] = {p0[9], p0[10], p0[8], p0[3], p0[4], p0[5]};
    // eps variation
    double p2[6] = {p0[9], p0[10], p0[2], p0[3], p0[4], 0.1};
    double p3[6] = {p1[0], p1[1], p1[2], p1[3], p1[4], 0.1};
    double p4[6] = {p0[9], p0[10], p0[2], p0[3], p0[4], 1.9};
    double p5[6] = {p1[0], p1[1], p1[2], p1[3], p1[4], 1.9};*/

    // double p3[6] = {p1[0], p1[1], p1[8], p1[3], p1[4], p0[5]};
    //
    // double *pa[n] = {p8, p2};
    // double *pa[n] = {p0arr, p1, p2, p3};
    //double *pa[n] = {p0arr};//, p1, p2, p3, p4, p5};
    //   double *pa[] = {p0.data(), p1, p2, p3};
    // double* pa[] = { p2, p4 };
    //  ceres::Problem *problem = new ceres::Problem[n]; // destructor segfaults
    ceres::Problem problem[n];
    for (unsigned i = 0; i < xdata.size(); i++)
        for (unsigned k = 0; k < n; k++) {
            ceres::CostFunction *cost_function =
                    new ceres::AutoDiffCostFunction<LossFunction, 2, 6>(
                            new LossFunction(xdata(i), ydata(i)));
            problem[k].AddResidualBlock(cost_function, nullptr, pa[k]);
        }
    // bounds
    for (unsigned i = 0; i < n; i++) {
        for (int j = 0; j < _parameter_count; j++) {
            problem[i].SetParameterLowerBound(pa[i], j, _parameter_mins[j]);
            problem[i].SetParameterUpperBound(pa[i], j, _parameter_maxs[j]);
        }
        // avoids auto-occlusion
        double xoff, yoff;
        if (abs(xdata.mean() - a._pose[0]) < 1)
            xoff = xdata.mean() / 2 + a._pose[0] / 2;
        else xoff = xdata.mean();
        if (abs(ydata.mean() - a._pose[1]) < 1)
            yoff = ydata.mean() / 2 + a._pose[1] / 2;
        else yoff = ydata.mean();
        if (a._pose(0) - xdata.mean() > 0) // robot on right of shape
            problem[i].SetParameterUpperBound(pa[i], 0, xoff); // x
        else
            problem[i].SetParameterLowerBound(pa[i], 0, xoff);
        if (a._pose(1) - ydata.mean() > 0)                               // robot above shape
            problem[i].SetParameterUpperBound(pa[i], 1, yoff); // y
        else
            problem[i].SetParameterLowerBound(pa[i], 1, yoff);
    }
    ceres::Solver::Options options;
    options.num_threads = _options.nthreads;
    options.minimizer_progress_to_stdout = false;
    options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY; //<- slower
    options.preconditioner_type = ceres::CLUSTER_JACOBI;       //<- also slower
    options.initial_trust_region_radius = 1e9;                 // important
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
        bool arrayIsFinite = true;
        for (int k = 0; k < 6; k++)
            if (!isfinite(pa[i][k])) {
                arrayIsFinite = false;
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
    e._p = ans;
    if (alreadyInitialized) ap_dop(e, a);
    else dop(e, a);
}

struct EllipseModel::CostFunction {
    CostFunction(double x_p,
                 double y_p,
                 double ap_d,
                 double an,
                 double sigma_d,
                 VectorXd ap_p,
                 MatrixXd ap_dop,
                 MatrixXd Q_a)
            : _ap_d(ap_d), _an(an), _x_p(x_p), _y_p(y_p), _sigma_d(sigma_d),
              _ap_p(ap_p), _ap_dop(ap_dop), Q(Q_a) {
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
        residual[1] = 50. * abs((d[0] - _ap_d) * (d[0] - _ap_d) / _sigma_d);
        // residual[1] *= 100000;
        //   parameter variation constraint
        Eigen::Map<const Eigen::Vector<T, 6>> pv(p);
        Eigen::Vector<T, 6> dp = pv - _ap_p;
        residual[2] = 30. * abs((dp.transpose() * _ap_dop.inverse() * Q.inverse() * dp)(0));
        // residual[2] *= 100000;
        //   area
        residual[3] = abs(pow((p[4]) * (p[3]), 0.5));
        // residual[3] = 0.;
        return true;
    }

private:
    double _ap_d, _an, _x_p, _y_p, _sigma_d;
    VectorXd _ap_p;
    MatrixXd _ap_dop, Q;
};

void EllipseModel::ap_ls(Entity &e, const Aggregate &a) const {
    const auto dan = a.get_rotated_measurement_mat();
    const auto xy = a.get_xy_mat();
    VectorXd d = dan.col(0);
    VectorXd an = dan.col(1);
    VectorXd x = xy.col(0);
    VectorXd y = xy.col(1);
    // TODO: GENERALIZE FOR N-DIM MEASUREMENTS
    ceres::Problem problem;
    // parameter count is the number of model params + number of measurements
    double n_p[6];
    vector<double> n_d(d.size());
    // fill initial values
    for (int i = 0; i < 6; i++)
        n_p[i] = e._p(i);
    for (int i = 0; i < d.size(); i++)
        n_d[i] = d[i];
    // fill problem
    for (unsigned i = 0; i < d.size(); i++) {
        // each parameter block optimizes the 6 parameters and one measurement
        // cout << endl << ap_p << endl << ap_dop << endl;
        ceres::CostFunction *cost_function =
                new ceres::AutoDiffCostFunction<CostFunction, 4, 6, 1>(
                        new CostFunction(a._pose(0), a._pose(1), d(i), an(i),
                                         _dop_sigma, e._p, e._E, _Q_a));
        problem.AddResidualBlock(cost_function, nullptr, n_p, &n_d[i]);
    }
    // bounds
    /*problem.SetParameterLowerBound(n_p, 3, 0.25); // a
    prblem.SetParameterLowerBound(n_p, 4, 0.25); // b
    problem.SetParameterLowerBound(n_p, 5, 0.3);  // eps <- very important
    problem.SetParameterUpperBound(n_p, 3, 30);   // a
    problem.SetParameterUpperBound(n_p, 4, 30);   // b
    problem.SetParameterUpperBound(n_p, 5, 1.7);*/  // eps
    for (int j = 2; j < _parameter_count - 1; j++) { // we skip x,y and eps
        problem.SetParameterLowerBound(n_p, j, _parameter_mins[j]);
        problem.SetParameterUpperBound(n_p, j, _parameter_maxs[j]);
    }
    problem.SetParameterUpperBound(n_p, 5, 1.9);
    problem.SetParameterLowerBound(n_p, 5, 0.1);
    // avoids auto-occlusion, only if enough points
    // if enough points we assume that we see enough of the ellipse
    // for an accurate mean
    // TODO: change to N
    if (x.size() > _options.init_npoints) {
        if (a._pose[0] - x.mean() > 0) // robot on right of shape
            problem.SetParameterUpperBound(n_p, 0, x.mean() / 2 + a._pose[0] / 2); // x
        else
            problem.SetParameterLowerBound(n_p, 0, x.mean() / 2 + a._pose[0] / 2);
        if (a._pose[1] - y.mean() > 0)                          // robot above shape
            problem.SetParameterUpperBound(n_p, 1, y.mean() / 2 + a._pose[1] / 2); // y
        else
            problem.SetParameterLowerBound(n_p, 1, y.mean() / 2 + a._pose[1] / 2);
    }
    // solve
    ceres::Solver::Options options;
    options.num_threads = _options.nthreads;
    options.minimizer_progress_to_stdout = false;
    // options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY; //<- slower
    // options.preconditioner_type = ceres::CLUSTER_JACOBI;       //<- also slower
    // options.initial_trust_region_radius = 1e8;                 // important
    options.max_num_iterations = 10;
    // Set stricter stopping criteria for more accurate results
    options.function_tolerance = 1e-14;    // Default: 1e-6
    options.gradient_tolerance = 1e-14;   // Default: 1e-10
    options.parameter_tolerance = 1e-14;   // Default: 1e-8

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // fixes
    VectorXd ans = Map<Vector<double, 6>>(n_p);
    //ans.conservativeResize(6 + d.size());
    Aggregate n_a(a);
    for (int i = 0; i < d.size(); i++)
        n_a._data_vector[i]._measurement(0) = n_d[i];
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
    e._p = ans(seq(0, 5));
    ap_dop(e, n_a);
}

bool EllipseModel::safety(Entity &e) const {
    if (e._p(5) < 0.1) {
        e._p(5) = 0.1;
        return true;
    }
    if (e._p(5) > 1.9) {
        e._p(5) = 1.9;
        return true;
    }
    if (e._p(3) < 0.25) {
        e._p(3) = 0.25;
        return true;
    }
    if (e._p(4) < 0.25) {
        e._p(4) = 0.25;
        return true;
    }
    if (e._p(2) > 2 * M_PI) {
        e._p(2) -= 2 * M_PI;
        return false;
    }
    if (e._p(2) < 0) {
        e._p(2) += 2 * M_PI;
        return false;
    }
    return false;
}

void EllipseModel::init_post(Entity &e, const Aggregate &a) const {}

bool EllipseModel::associate(const Entity &e, const Data &d) const {
    return true;
}

void EllipseModel::augment_post(Entity &e, const Data &d) const {}


double EllipseModel::fss(const VectorXd &p, const VectorXd &pose, const VectorXd &rotated_measurement) const {
    return ft<double>(
            p(0), p(1), p(2), p(3), p(4), p(5),
            pose(0), pose(1), rotated_measurement(0), rotated_measurement(1),
            0, 0);
}

void EllipseModel::decimate(Entity &e) const {
    if (!e._a) return;
    // we will first get all angles for the points with respect to the center of the ellipse
    Aggregate newAggregate;
    auto pxy = e._a->get_xy_mat();
    VectorXd angles = atan2((pxy.col(1).array() - e._p(1)),
                            (pxy.col(0).array() - e._p(0))).array() + M_PI;
    int lastInRangeIdx = 0;
    Aggregate slices[SLICE_NUM];
    //double inc = M_PI / 100;
    for (int i = 0; i < angles.size(); i++) {
        int index = static_cast<int>(floor(angles(i) / (2 * M_PI / SLICE_NUM))) % SLICE_NUM;
        slices[index].push_back(e._a->_data_vector[i]);
    }
    //for(int i=0;i<SLICE_NUM;i++)
    //    cout<<"SLICE "<<i<<" has "<< slices[i]._data_vector.size()<<" datas"<<endl;

    // random sampling
    // samples per slice
    random_device rd;
    mt19937 gen(rd());
    for (auto &slice: slices) {
        uniform_int_distribution<> distrib(0, slice._data_vector.size() - 1);
        int N = min(3, (int) slice._data_vector.size()); // samples per slice
        // Generate a vector of indices
        vector<int> indices(slice._data_vector.size());
        iota(indices.begin(), indices.end(), 0);
        shuffle(indices.begin(), indices.end(), gen);
        for (int k = 0; k < N; k++)
            newAggregate.push_back(slice._data_vector[k]);
    }
    /*for (double b = 0; b < 2 * M_PI; b += inc) {
        while (angles(lastInRangeIdx) < b + inc &&
               angles(lastInRangeIdx) >= b &&
               lastInRangeIdx < angles.size() - 1)
            lastInRangeIdx++;
        newAggregate.push_back(e._a->_data_vector[lastInRangeIdx]);
    }*/
    e._a->clear();
    e._a->push_back(newAggregate);
}

template<typename T>
T EllipseModel::ft(T xc, T yc, T th, T a, T b, T e, T xp, T yp, T d, T an, T mud, T muan) const {
    T x = xp + (d + mud) * cos(an + muan);
    T y = yp + (d + mud) * sin(an + muan);
    T f1 = ((x - xc) * cos(th) + (y - yc) * sin(th)) / a;
    T f2 = ((x - xc) * sin(th) - (y - yc) * cos(th)) / b;
    T buf = pow(pow(f1, 2.), 1. / e) + pow(pow(f2, 2.), 1. / e);
    return sqrt(a * b) * (pow(buf, e) - 1.);
}

