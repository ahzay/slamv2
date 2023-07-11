//
// Created by user on 7/10/23.
//

#include "tools.h"
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
