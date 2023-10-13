#ifndef PREPROC_HPP
#define PREPROC_HPP

#include "visualizer.hpp"
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <iostream>

namespace fsm = boost::filesystem;
using namespace Eigen;

//
class Scan {
    // Matrix<double, Dynamic, 2> pts, mes;
public:
    Matrix<double, Dynamic, 4> data; // x, y, d, an
    Vector<double, 2> loc;
    Vector<double, 2> odomloc;
    double ori;

    void write_scan(int n) {
        ofstream f("scan_" + to_string(n) + ".csv");
        // first line: posx,posy,ori,rows
        f << loc(0) << "," << loc(1) << "," << ori << "," << data.rows() << endl;
        // subsequent lines: x,y,d,an
        for (int i = 0; i < data.rows(); i++)
            f << data(i, 0) << "," << data(i, 1) << "," << data(i, 2) << ","
              << data(i, 3) << endl;
        f.close();
    }

    void read_scan(const string &dir, int n, int z,
                   float mult) { // z: to remove from each side
        char cbuf;
        double dbuf;
        string filename = dir + "/combined_" + to_string(n) + ".csv";
        ifstream f(filename);
        if (!f.is_open()) {
            throw std::runtime_error("Error opening the scan file!");
        }
        int rows = 0;
        // first line
        f >> loc(0) >> cbuf >> loc(1) >> cbuf >> ori;
        // mult
        loc *= mult;
        // odom x y z q0 q1 q2 q3
        f >> odomloc(0) >> odomloc(1) >> dbuf >> dbuf >> dbuf >> dbuf >> dbuf;
        // ,rows
        f >> cbuf >> rows;
        cout << "loc: " << loc.transpose() << " ori: " << ori << endl;
        data.conservativeResize(0, NoChange);
        // get data
        for (int i = 0; i < rows; i++) {
            long j = data.rows();
            Vector<double, 4> bvec;
            f >> bvec(0) >> cbuf >> bvec(1) >> cbuf >> bvec(2) >> cbuf >> bvec(3);
            if (bvec(2) > 0 && bvec(2) < 5) { // restriction on min/max distance
                data.conservativeResize(data.rows() + 1, NoChange);
                data.row(j) = bvec;
            }
        }
        // mult and error
        for (int i = 0; i < data.rows(); i++) {
            double error = (-0.01 + (rand() / (RAND_MAX / 0.02)));
            // cout << "error: " << error << endl;
            data.col(2)(i) += error;
        }
        data.col(2) *= mult;
        // data.col(2).array() *= 3; // fixing the scan ...
        data = data(seq(z, data.rows() - z - 1), all);
        // rows = data.rows();
        //  adjust data to emulate null rotation
        data.col(3).array() += ori;
        data.col(0) = data.col(2).array() * data.col(3).array().cos() + loc(0);
        data.col(1) = data.col(2).array() * data.col(3).array().sin() + loc(1);
    }
};

// angle constraint 0, 2pi
double ancnstr(double x) {
    x = fmod(x, 2 * M_PI);
    if (x < 0)
        x += 2 * M_PI;
    return x;
}

bool compare_angle(const VectorXd &lhs, const VectorXd &rhs) {
    return lhs(3) < rhs(3);
}

#endif // PREPROC_HPP
