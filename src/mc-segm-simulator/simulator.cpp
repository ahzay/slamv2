//
// Created by user on 10/11/23.
//

#include "simulator.h"
#include <fstream>
#include <utility>
#include <cmath>
#include <omp.h>

Simulator::Simulator(CmdLineOptions options) : o(std::move(options)) {
    generator.seed(std::random_device()());
}

double Simulator::get_rand(double min, double max) {
    uniform_real_distribution<double> distribution(min, max);
    return distribution(generator);
}

Vector<double, 6> Simulator::gen_params() {
    Vector<double, 6> p;
    p(0) = get_rand(o.xmin, o.xmax);
    p(1) = get_rand(o.ymin, o.ymax);
    p(2) = get_rand(-M_PI, M_PI);
    p(3) = get_rand(0.5, (o.xmax - o.xmin) / 4);
    p(4) = get_rand(0.5, (o.ymax - o.ymin) / 4);
    p(5) = get_rand(0.1, 1.9);
    return p;
    //ofstream of("ground_truth.txt");
    //of<<p.transpose();
}

bool Simulator::detect_overlap(Vector<double, 6> p1, Vector<double, 6> p2) {
    double dx = p1[0] - p2[0];
    double dy = p1[1] - p2[1];
    double distance = std::sqrt(dx * dx + dy * dy);
    double radius1 = std::max(p1[2], p1[3]);
    double radius2 = std::max(p2[2], p2[3]);
    return distance < (radius1 + radius2);
}


void Simulator::gen_all_params() {
    ps.push_back(gen_params());
    while (ps.size() < o.nell) {
        auto testp = gen_params();
        if (!any_of(ps.begin(), ps.end(), [&testp, this](const Vector<double, 6> &p) {
            return detect_overlap(testp, p);
        }))
            ps.push_back(testp);
    }
}

void Simulator::gen_location() {
    do {
        l(0) = get_rand(o.xmin, o.xmax);
        l(1) = get_rand(o.ymin, o.ymax);
    } while (std::any_of(ps.begin(), ps.end(), [this](const Vector<double, 6> &p) {
        return sqrt(pow(l(0) - p(0), 2) + pow(l(1) - p(1), 2)) <= max(p(3), p(4));
    }));
}

void Simulator::gen_all() {
    gen_all_params();
    gen_location();
    ofstream of("ground_truth.txt");
    if (of.is_open()) {
        of << ps.size() << endl;  // Write the number of ellipses to the first line
        for (const auto &p: ps)
            of << p.transpose() << endl;  // New line for each ellipse
    }
}

void Simulator::simulate() {
    gen_all();
    const double dmax = 20;//(l({0, 1}) - p({0, 1})).norm() + max(p(3), p(4));
    const int max_size = static_cast<int>(2 * M_PI / ((o.an / o.npts) * M_PI / 180.0));
    ms.resize(max_size);  // Pre-allocate space
    MWI buf = {Vector2d(0, 0), -1};
    ms.assign(max_size, buf);
#pragma omp parallel for default(none) shared(ms, ps, l, o, max_size, dmax)
    for (int i = 0; i < max_size; ++i) {
        double an = i * (o.an / o.npts) * M_PI / 180.0;
        double e = MAXFLOAT;
        double var_e = 1000000;
        double d = 0;
        int min_p_index = -1;
        while (var_e > 1e-24 && e > 0 && d <= dmax) {
            double x = l(0) + d * cos(an);
            double y = l(1) + d * sin(an);
            std::vector<double> e_values;
            std::vector<Vector<double, 6>>::iterator min_p_it = ps.begin(); // Iterator to the p that gives min e
            double min_e = std::numeric_limits<double>::max(); // Initialize with max double value
            for (auto it = ps.begin(); it != ps.end(); ++it) {
                double e_val = calc_e(*it, x, y);
                e_values.push_back(e_val);
                if (e_val < min_e) {
                    min_e = e_val;
                    min_p_it = it;
                }
            }
            double new_e = min_e; // The minimum e value
            var_e = abs(new_e - e);
            e = new_e;
            d += std::min(0.01, e);
            min_p_index = static_cast<int>(std::distance(ps.begin(), min_p_it));
        }
        if (d <= dmax && e <= 0) {
            ms[i] = {Vector2d(d, an), min_p_index};  // Direct assignment
        }
    }
    // remove empty
    ms.erase(std::remove_if(ms.begin(), ms.end(), [](const auto &m) {
        return m.m(0) == 0.0 && m.m(1) == 0.0;
    }), ms.end());
    // writing phase
    ofstream of(o.data_folder + "/combined_0.csv");
    if (!of.is_open()) {
        throw std::runtime_error("Failed to open the file, check that the folder exists");
    }
    of << l(0) << "," << l(1) << "," << "0.0" << endl;
    of << "0.0 0.0 0.0 0.0 0.0 0.0 0.0" << endl;
    of << "," << ms.size() << endl;
    for (const auto &m: ms)
        of << "0.0,0.0," << m.m(0) << "," << m.m(1) << endl;
    of.close();
    of.open("ground_segm_truth.txt");
    if (of.is_open()) {
        of << ms.size() << endl;
        for (const auto &m: ms)
            of << m.m.transpose() << " " << m.i << endl;
    }
}

double Simulator::calc_e(Vector<double, 6> p, double x, double y) {
    double f1 = ((x - p(0)) * cos(p(2)) + (y - p(1)) * sin(p(2))) / p(3);
    double f2 = ((x - p(0)) * sin(p(2)) - (y - p(1)) * cos(p(2))) / p(4);
    double buf = pow(pow(f1, 2.), 1. / p(5)) + pow(pow(f2, 2.), 1. / p(5));
    double e = sqrt(p(3) * p(4)) * (pow(buf, p(5)) - 1.);
    return e;
}





