//
// Created by user on 10/11/23.
//

#include "simulator.h"
#include <fstream>
#include <utility>
#include <cmath>

Simulator::Simulator(CmdLineOptions options) : o(std::move(options)) {
    generator.seed(std::random_device()());
}

double Simulator::get_rand(double min, double max) {
    uniform_real_distribution<double> distribution(min, max);
    return distribution(generator);
}

void Simulator::gen_params() {
    p(0) = get_rand(o.xmin, o.xmax);
    p(1) = get_rand(o.ymin, o.ymax);
    p(2) = get_rand(-M_PI, M_PI);
    p(3) = get_rand(0.5, (o.xmax - o.xmin) / 4);
    p(4) = get_rand(0.5, (o.ymax - o.ymin) / 4);
    p(5) = get_rand(0.1, 1.9);
    ofstream of("ground_truth.txt");
    of<<p.transpose();
}

void Simulator::gen_location() {
    // exclusion zone circle radius
    double r = max(p(3), p(4));
    do {
        l(0) = get_rand(o.xmin, o.xmax);
        l(1) = get_rand(o.ymin, o.ymax);
    } while (sqrt(pow(l(0) - p(0), 2) + pow(l(1) - p(1), 2)) <= r);
}

void Simulator::simulate() {
    gen_params();
    gen_location();
    double heading = atan2(p(1) - l(1), p(0) - l(0));
    double dmax = (l({0, 1}) - p({0, 1})).norm() + max(p(3), p(4));
    for (double an = heading - M_PI / 2; an < heading + M_PI / 2; an += (o.an / o.npts) * M_PI / 180.0) {
        double e = MAXFLOAT;
        double var_e = 1000000;
        double d = 0;
        //cout << "FOR AN: " << an << endl;
        while (var_e > 1e-24 && e >= 0 && d <= dmax) {
            // calculate error
            double x = l(0) + d * cos(an);
            double y = l(1) + d * sin(an);
            double f1 = ((x - p(0)) * cos(p(2)) + (y - p(1)) * sin(p(2))) / p(3);
            double f2 = ((x - p(0)) * sin(p(2)) - (y - p(1)) * cos(p(2))) / p(4);
            double buf = pow(pow(f1, 2.), 1. / p(5)) + pow(pow(f2, 2.), 1. / p(5));
            double new_e = sqrt(p(3) * p(4)) * (pow(buf, p(5)) - 1.);
            var_e = abs(new_e - e);
            e = new_e;
            d += min(0.01, e);
            //cout << d << endl;
        }
        if (d <= dmax)
            ms.push_back(Vector2d(d, an));
    }
    // writing phase
    ofstream of(o.data_folder + "/combined_0.csv");
    if (!of.is_open()) {
        throw std::runtime_error("Failed to open the file, check that the folder exists");
    }
    of << l(0) << "," << l(1) << "," << "0.0" << endl;
    of << "0.0 0.0 0.0 0.0 0.0 0.0 0.0" << endl;
    of << "," << ms.size() << endl;
    for (const auto &m: ms)
        of << "0.0,0.0," << m(0) << "," << m(1) << endl;

}


