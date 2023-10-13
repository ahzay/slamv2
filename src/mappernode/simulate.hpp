#ifndef SIMULATE_HPP
#define SIMULATE_HPP

#include "vecext.hpp"
#include <omp.h>
#include <random>
#include <tuple>

struct Object {
    vecext<double> p; // parameters
    unsigned t;       // 0 for ellipse (6 params), 1 for line (4 params)
};

template<typename T>
class Simulator {
public:
    tuple<vecext<T>, vecext<T>, vecext<T>, vecext<T>>
    simulate(vecext<T> p, vecext<T> loc, T theta, T rnd) {
        random_device rd;  // used to obtain a seed for the random number engine
        mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        uniform_real_distribution<> dist(0.0, 1.0);
        // optimizing t to only look in the direction of the object
        double direction = atan2(p[1] - loc[1], p[0] - loc[0]);
        vecext<T> t = linspace(direction - M_PI / 2, direction + M_PI / 2,
                               size_t(2 * M_PI / theta) + 1);
        vecext<T> x, y, d, th;
#pragma omp critical
#pragma omp parallel for ordered shared(x, y)
        for (T angle: t) {
            vecext<T> X = linspace(0.0, 40.0, 100000) * cos(angle) + loc[0];
            vecext<T> Y = linspace(0.0, 40.0, 100000) * sin(angle) + loc[1];
            vecext<T> res = pow(pow(f1(p, X, Y), 2.0), (1 / p[5])) +
                            pow(pow(f2(p, X, Y), 2.0), (1 / p[5])) - 1.0;
            int min_idx = res.minimas_idx(0);
            T min = 10;
            if (min_idx != -1)
                min = abs(res.at(min_idx));
#pragma omp ordered
            if (min < 0.01) {
                // TODO: fix this so x and y correspond to d
                d.push_back(linspace(0.0, 40.0, 100000)[min_idx] + dist(gen) * rnd -
                            rnd / 2);
                th.push_back(angle);
                x.push_back(d.back() * cos(th.back()) + loc[0]);
                y.push_back(d.back() * sin(th.back()) + loc[1]);
            }
        }
        // for simulation glitches
        double verif = (pow(pow(f1(p, x, y), 2.0), (1 / p[5])) +
                        pow(pow(f2(p, x, y), 2.0), (1 / p[5])) - 1.0)
                               .sum() /
                       x.size();
        if (abs(verif) < 0.1)
            return {x, y, d, th};
        else
            exit(1);
    }

    tuple<vecext<T>, vecext<T>, vecext<T>, vecext<T>>
    simulate(tuple<vecext<T>, vecext<T>, vecext<T>, vecext<T>> data, vecext<T> p,
             vecext<T> loc, T theta, T rnd) {
        tuple<vecext<T>, vecext<T>, vecext<T>, vecext<T>> datanew =
                simulate(p, loc, theta, rnd);
        get<0>(data).insert(get<0>(data).end(), get<0>(datanew).begin(),
                            get<0>(datanew).end());
        get<1>(data).insert(get<1>(data).end(), get<1>(datanew).begin(),
                            get<1>(datanew).end());
        get<2>(data).insert(get<2>(data).end(), get<2>(datanew).begin(),
                            get<2>(datanew).end());
        get<3>(data).insert(get<3>(data).end(), get<3>(datanew).begin(),
                            get<3>(datanew).end());
        return data;
    }

    tuple<vecext<T>, vecext<T>, vecext<T>, vecext<T>>
    simulate2(vector<vecext<T>> ps, vecext<T> loc, T theta, T rnd) {

        random_device rd;  // used to obtain a seed for the random number engine
        mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        uniform_real_distribution<> dist(0.0, 1.0);
        // optimizing t to only look in the direction of the object
        // double direction = atan2(p[1] - loc[1], p[0] - loc[0]);
        vecext<T> t = linspace(0.0, 2 * M_PI, size_t(2 * M_PI / theta) + 1);
        vecext<T> x, y, d, th;

        //#pragma omp critical
        //#pragma omp parallel for ordered shared(x, y)
        for (T angle: t) {
            for (auto &p: ps) {
                vecext<T> X = linspace(0.0, 40.0, 100000) * cos(angle) + loc[0];
                vecext<T> Y = linspace(0.0, 40.0, 100000) * sin(angle) + loc[1];
                vecext<T> res = pow(pow(f1(p, X, Y), 2.0), (1 / p[5])) +
                                pow(pow(f2(p, X, Y), 2.0), (1 / p[5])) - 1.0;
                int min_idx = res.minimas_idx(0);
                T min = 10;
                if (min_idx != -1)
                    min = abs(res.at(min_idx));
                //#pragma omp ordered
                if (min < 0.001) {
                    d.push_back(linspace(0.0, 40.0, 100000)[min_idx] + dist(gen) * rnd -
                                rnd / 2);
                    th.push_back(angle);
                    x.push_back(d.back() * cos(th.back()) + loc[0]);
                    y.push_back(d.back() * sin(th.back()) + loc[1]);
                    goto next;
                }
            }
            next:;
        }
        return {x, y, d, th};
    }

    tuple<vecext<double>, vecext<double>, vecext<double>, vecext<double>>
    simulate3(vector<Object> os, vecext<double> loc, double theta, double rnd) {

        random_device rd;  // used to obtain a seed for the random number engine
        mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        uniform_real_distribution<> dist(0.0, 1.0);
        // optimizing t to only look in the direction of the object
        // double direction = atan2(p[1] - loc[1], p[0] - loc[0]);
        vecext<T> t = linspace(0.0, 2 * M_PI, size_t(2 * M_PI / theta) + 1);
        vecext<T> x, y, d, th;
        double l = 10;
        for (T angle: t) {
            for (auto &o: os) {
                vecext<T> X = linspace(0.0, l, 100000) * cos(angle) + loc[0];
                vecext<T> Y = linspace(0.0, l, 100000) * sin(angle) + loc[1];

                if (o.t == 0) {
                    vecext<T> res = pow(pow(f1(o.p, X, Y), 2.0), (1 / o.p[5])) +
                                    pow(pow(f2(o.p, X, Y), 2.0), (1 / o.p[5])) - 1.0;
                    int min_idx = res.minimas_idx(0);
                    T min_num = 10;
                    if (min_idx != -1)
                        min_num = abs(res.at(min_idx));
                    if (min_num < 0.001) {
                        d.push_back(linspace(0.0, l, 100000)[min_idx] + dist(gen) * rnd -
                                    rnd / 2);
                        th.push_back(angle);
                        x.push_back(d.back() * cos(th.back()) + loc[0]);
                        y.push_back(d.back() * sin(th.back()) + loc[1]);
                        goto next;
                    }
                }
                if (o.t == 1) {
                    vecext<T> res = // p[1]=r,p[0]=a
                            (1.0 / o.p[1]) * (X * cos(o.p[0]) + Y * sin(o.p[0])) - 1.0;
                    int min_idx = res.minimas_idx(0);
                    T min_num = 10;
                    if (min_idx != -1)
                        min_num = abs(res.at(min_idx));
                    if (min_num < 0.001) {
                        double xtest =
                                linspace(0.0, l, 100000)[min_idx] * cos(angle) + loc[0];
                        double ytest =
                                linspace(0.0, l, 100000)[min_idx] * sin(angle) + loc[1];
                        double antest = atan2(ytest, xtest);
                        if (antest >= min(o.p[2], o.p[3]) &&
                            antest <= max(o.p[2], o.p[3])) {

                            d.push_back(linspace(0.0, l, 100000)[min_idx] + dist(gen) * rnd -
                                        rnd / 2);
                            th.push_back(angle);
                            x.push_back(d.back() * cos(th.back()) + loc[0]);
                            y.push_back(d.back() * sin(th.back()) + loc[1]);
                            goto next;
                        }
                    }
                }
            }
            next:;
        }
        return {x, y, d, th};
    }

private:
    vecext<T> f1(vecext<T> p, vecext<T> x, vecext<T> y) {
        return ((x - p[0]) * cos(p[2]) + (y - p[1]) * sin(p[2])) / p[3];
    }

    vecext<T> f2(vecext<T> p, vecext<T> x, vecext<T> y) {
        return ((x - p[0]) * sin(p[2]) - (y - p[1]) * cos(p[2])) / p[4];
    }

    vecext<T> fs(vecext<T> p, vecext<T> x, vecext<T> y) {
        vecext<T> j1 = pow(pow(f1(p, x, y), 2.0), (1 / p[5]));
        vecext<T> j2 = pow(pow(f2(p, x, y), 2.0), (1 / p[5]));
        return pow((j1 + j2), p[5]) - 1;
    }
};

vector<tuple<vecext<double>, vecext<double>, vecext<double>, vecext<double>>>
full_sim(vector<vecext<double>> ps, vector<vecext<double>> ls, double theta,
         double rnd) {
    // sim object
    Simulator<double> sim;
    // declare result variable
    vector<tuple<vecext<double>, vecext<double>, vecext<double>, vecext<double>>>
            res;
    // iterate over locations
    tuple<vecext<double>, vecext<double>, vecext<double>, vecext<double>> data;
    //#pragma omp critical
    //#pragma omp parallel for ordered shared(res)
    for (auto &l: ls) {
        // iterate over objects
        // auto data = sim.simulate(ps[0], l, theta, rnd);
        // for (unsigned i = 1; i < ps.size(); i++) {
        //  data = sim.simulate(data, ps[i], l, theta, rnd);
        //}
        data = sim.simulate2(ps, l, theta, rnd);
        //#pragma omp ordered
        { res.push_back(data); }
    }
    return res;
}

vector<tuple<vecext<double>, vecext<double>, vecext<double>, vecext<double>>>
full_sim2(vector<Object> os, vector<vecext<double>> ls, double theta,
          double rnd) {
    // sim object
    Simulator<double> sim;
    // declare result variable
    vector<tuple<vecext<double>, vecext<double>, vecext<double>, vecext<double>>>
            res;
    // iterate over locations
    tuple<vecext<double>, vecext<double>, vecext<double>, vecext<double>> data;
    //#pragma omp critical
    //#pragma omp parallel for ordered shared(res)
    for (auto &l: ls) {
        // iterate over objects
        // auto data = sim.simulate(ps[0], l, theta, rnd);
        // for (unsigned i = 1; i < ps.size(); i++) {
        //  data = sim.simulate(data, ps[i], l, theta, rnd);
        //}
        data = sim.simulate3(os, l, theta, rnd);
        //#pragma omp ordered
        { res.push_back(data); }
    }
    return res;
}

/*
class Simulator2 {
public:
  tuple<vecext<double>, vecext<double>> simulate(double *p, double *loc,
                                                 double theta, double eps) {
    Eigen::ArrayXd t =
        Eigen::ArrayXd::LinSpaced(size_t(2 * M_PI / theta) + 1, 0, 2 * M_PI);
    vecext<double> x, y;
    for (double a : t) {
      Eigen::ArrayXd X =
          Eigen::ArrayXd::LinSpaced(100000, 0, 30) * cos(a) + loc[0];
      Eigen::ArrayXd Y =
          Eigen::ArrayXd::LinSpaced(100000, 0, 30) * sin(a) + loc[1];
      Eigen::ArrayXd res = f1(p, X, Y).pow(2).pow(1 / p[5]) +
                           f2(p, X, Y).pow(2).pow(1 / p[5]) - 1;
      // finding min_idx
      int min_idx = -1;
      for (int i = 1; i < res.size() - 1; i++) {
        double prev = abs(res(i - 1));
        double curr = abs(res(i));
        double next = abs(res(i + 1));
        if (prev > curr && curr < next)
          min_idx = i;
      }
      // answers
      double min = 10;
      if (min_idx != -1)
        min = abs(res(min_idx));
      if (min < 0.01) {
        x.push_back(X(min_idx));
        y.push_back(Y(min_idx));
      }
    }
    return {x, y};
  }

private:
  Eigen::ArrayXd f1(double *p, Eigen::ArrayXd x, Eigen::ArrayXd y) {
    return ((x - p[0]) * cos(p[2]) + (y - p[1]) * sin(p[2])) / p[3];
  }
  Eigen::ArrayXd f2(double *p, Eigen::ArrayXd x, Eigen::ArrayXd y) {
    return ((x - p[0]) * sin(p[2]) - (y - p[1]) * cos(p[2])) / p[4];
  }
};
*/
#endif // SIMULATE_HPP
