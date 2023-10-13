//
// Created by user on 7/13/23.
//

#ifndef MAPPERNODE_LINEMODEL_H
#define MAPPERNODE_LINEMODEL_H

#include "model.h"
#include <autodiff/reverse/var.hpp>
#include <ceres/autodiff_cost_function.h>
#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>

namespace ad = autodiff;

class LineModel : public Model {
public:
    explicit LineModel(const string &file, const CmdLineOptions &options);

    Matrix<double, 1, -1> dfs(const Entity &e, const Data &d) const override;

    VectorXd init(const Aggregate &a) const;

    void ls(Entity &e, const Aggregate &a) const override;

    void ap_ls(Entity &e, const Aggregate &a) const override;

    bool safety(Entity &e) const override;

    void init_post(Entity &e, const Aggregate &a) const override;

    bool associate(const Entity &e, const Data &d) const override;

    void augment_post(Entity &e, const Data &d) const override;

    template<typename T>
    T ft(T r, T a, T xp, T yp, T d, T an, T mud, T muan) const;

    double fss(const VectorXd &p, const VectorXd &pose, const VectorXd &rotated_measurement) const override;

    struct CostFunction;
};


#endif //MAPPERNODE_LINEMODEL_H
