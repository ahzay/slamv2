//
// Created by user on 7/10/23.
//

#include "ellipsemodel.h"
#include "entity.h"
EllipseModel::EllipseModel(const string &file) : Model(file) {}

VectorXd EllipseModel::fs(const Entity &e, const Aggregate &a) {
    const auto data = a.get_mat();
    VectorXd ans(data.rows());
    for (int i = 0; i < data.rows(); i++)
        ans(i) = fss<double>(e._p, a._pos, data(i, all));
    return ans;
}

MatrixXd EllipseModel::dfs(const Entity &e, const Aggregate &a) {
    const auto data = a.get_mat();
    ad::var xc(e._p(0)), yc(e._p(1)), th(e._p(2)), aa(e._p(3)), b(e._p(4)), ee(e._p(5)), xp(a._pos(0)),
            yp(a._pos(1)), d(0.), an(0.), mud(0.), muan(0.);
    ad::var f = ft(xc, yc, th, aa, b, ee, xp, yp, d, an, mud, muan);

    MatrixXd ans(data.rows(), 12);
    for (unsigned i = 0; i < data.rows(); i++) {
        d.update(data.row(i)(2));
        an.update(data.row(i)(3));
        f.update();
        auto [dxc, dyc, dth, da, db, de, dxp, dyp, dd, dan, dmud, dmuan] =
                ad::derivatives(f,
                                ad::wrt(xc, yc, th, aa, b, ee, xp, yp, d, an, mud, muan));
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
    }
    return ans;
}

template<typename T>
T EllipseModel::fss(const VectorXd &p, const VectorXd &pos, const VectorXd &data) {
    return ft<T>(
            p(0), p(1), p(2), p(3), p(4), p(5), pos(0), pos(1), data(2), data(3), 0, 0);
}

template<typename T>
T EllipseModel::ft(T xc, T yc, T th, T a, T b, T e, T xp, T yp, T d, T an, T mud, T muan) {
    T x = xp + (d + mud) * cos(an + muan);
    T y = yp + (d + mud) * sin(an + muan);
    T f1 = ((x - xc) * cos(th) + (y - yc) * sin(th)) / a;
    T f2 = ((x - xc) * sin(th) - (y - yc) * cos(th)) / b;
    T buf = pow(pow(f1, 2.), 1. / e) + pow(pow(f2, 2.), 1. / e);
    return sqrt(a * b) * (pow(buf, e) - 1.);
}
