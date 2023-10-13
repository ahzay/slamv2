//
// Created by user on 7/10/23.
//

#include "tools.h"


int
sgn(double val) {
    return (double(0) < val) - (val < double(0));
}

double angle_diff(double angle1, double angle2) {
    const double pi = M_PI;
    const double twoPi = 2.0 * pi;
    double diff = fmod(abs(angle1 - angle2), twoPi);
    if (diff > pi) {
        diff = twoPi - diff;
    }
    return diff;
}


VectorXd
atan2(VectorXd y, VectorXd x) {
    VectorXd res(y.size());
    for (int i = 0; i < y.size(); i++) {
        res(i) = atan2(y(i), x(i));
    }
    return res;
}

// TODO: eliminates entity instead of stopping entire program
bool is_psd(const MatrixXd &A) {
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

double vector_cov(VectorXd x, VectorXd y) {
    auto xm = x.array() - x.mean();
    auto ym = y.array() - y.mean();
    auto m = xm.array() * ym.array();
    return m.sum() / (m.size() - 1);
}