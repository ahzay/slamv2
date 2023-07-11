//
// Created by user on 7/10/23.
//

#ifndef MAPPERNODE_ELLIPSEMODEL_H
#define MAPPERNODE_ELLIPSEMODEL_H

#include "model.h"
#include <autodiff/reverse/var.hpp>
#include <ceres/autodiff_cost_function.h>
#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>
namespace ad=autodiff;
class EllipseModel : public Model {
public:
    EllipseModel(const string &file);

    VectorXd fs(const Entity &e, const Aggregate &a) override;

    MatrixXd dfs(const Entity &e, const Aggregate &a) override;

    VectorXd ls(const Aggregate &a) override;

    VectorXd ls(Entity &e, const Aggregate &a) override;

    bool safety(Entity &e) override;

    void init_post(Entity &e, const Aggregate &a) override;

    bool associate(const Entity &e, const Data &d) override;

    void augment_post(Entity &e, const Data &d) override;

private:
    template<typename T>
    T ft(T xc, T yc, T th, T a, T b, T e, T xp, T yp, T d, T an, T mud, T muan) ;
    template<typename T>
    T fss(const VectorXd& p, const VectorXd &pos, const VectorXd &data);
};


#endif //MAPPERNODE_ELLIPSEMODEL_H
