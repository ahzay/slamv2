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

namespace ad = autodiff;

class EllipseModel : public Model {
public:
    explicit EllipseModel(const string &file);

    double fs(const Entity &e, const Data &d) const override;

    VectorXd dfs(const Entity &e, const Data &d) const override;

    VectorXd init(const Aggregate &a) const;

    VectorXd ls(const Aggregate &a) const override;

    VectorXd ls(Entity &e, const Aggregate &a) const override;

    bool safety(Entity &e) const override;

    void init_post(Entity &e, const Aggregate &a) const override;

    bool associate(const Entity &e, const Data &d) const override;

    void augment_post(Entity &e, const Data &d) const override;

private: // gross wrappers
    template<typename T>
    T ft(T xc, T yc, T th, T a, T b, T e, T xp, T yp, T d, T an, T mud, T muan) const;

    template<typename T>
    T fss(const VectorXd &p, const VectorXd &pose, const VectorXd &rotated_measurement) const;

    struct LossFunction;
    struct CostFunction;
};


#endif //MAPPERNODE_ELLIPSEMODEL_H
