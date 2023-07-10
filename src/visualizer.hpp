#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP

#include "simulate.hpp"
#include <Eigen/Core>
#include <fstream>

class Visualizer {
public:
    Visualizer(int cnt);

    void add_points(Eigen::Matrix<double, Eigen::Dynamic, 2> pts, const string& color);

    void add_ellipse(Eigen::VectorXd p);

    void add_segment(Eigen::VectorXd p, Eigen::VectorXd t);

    void save();

    ofstream of;
};

Visualizer::Visualizer(int cnt) {
    of.open(to_string(cnt) + ".gpt");
    of << "set term pngcairo size 1000,2000 enhanced font \"Times,12\" "
       << endl
       //<< "set term pngcairo size 1000,1000 enhanced font \"Times,12\" " << endl
       << "set size ratio -1 "
       << endl
       //<< "set contour" << endl
       << "set multiplot"
       << endl
       //<< "set cntrparam levels discrete 0 " << endl
       //<< "set view map" << endl
       //<< "unset surface " << endl
       << "set nokey "
       << endl
       //<< "set dgrid3d 1000,1000 qnorm 2"
       << endl
       << "set samples 2000" << endl
       << "set isosamples 10000,10000 " << endl
       << "set xrange [-12:12] " << endl
       << "set yrange [-21:21] "
       //<< endl
       //<< "set xrange [-5:5] " << endl
       //<< "set yrange [-15:-5] "
       << endl
       //<< "set xrange [-8:1] " << endl
       //<< "set yrange [-21:-12] " << endl
       << "set ytics 0.5" << endl
       << "set xtics 0.5" << endl
       << "set grid" << endl
       << "set grid mytics" << endl
       << "set grid mxtics"
       << endl
       //<< "set lmargin 0" << endl
       //<< "set bmargin 0" << endl
       //<< "set tmargin 0" << endl
       //<< "set rmargin 0" << endl
       << "se_x(t) = xc + a * cos(th) * sgn(cos(t)) * abs(cos(t))**(e) - b "
          "* sin(th) * sgn(sin(t)) * abs(sin(t))**(e)"
       << endl
       << "se_y(t) = yc + a * sin(th) * sgn(cos(t)) * abs(cos(t))**(e) + b "
          "* cos(th) * sgn(sin(t)) * abs(sin(t))**(e)"
       << endl
       << "set parametric" << endl
       << "set trange [0.001:2*pi-0.001]" << endl;
}

void Visualizer::add_points(Eigen::Matrix<double, Eigen::Dynamic, 2> pts,
                            const string& color) {
    // of << "unset border\n unset xtics\n unset ytics\n";
    of << R"(plot "-" w p ls 7 lw 0.05 lc rgb ")" << color << "\"" << endl;
    for (int i = 0; i < pts.rows(); i++)
        of << pts(i, 0) << " " << pts(i, 1) << endl;
    of << "e" << endl;
    // of << "set border\n set xtics\n set ytics\n";
}

void Visualizer::add_ellipse(Eigen::VectorXd p) {
    of << "xc=" << p(0) << endl
       << "yc=" << p(1) << endl
       << "th=" << p(2) << endl
       << "a=" << p(3) << endl
       << "b=" << p(4) << endl
       << "e=" << p(5) << endl
       << "plot se_x(t), se_y(t) with l ls 7 lw 3 lc rgb \"green\" "
          "notitle"
       << endl;
}

void Visualizer::add_segment(Eigen::VectorXd p, Eigen::VectorXd t) {
    double z_0 = p(0) / cos(abs(t(0) - p(1)));
    double z_1 = p(0) / cos(abs(t(1) - p(1)));
    double x_0 = z_0 * cos(t(0));
    double y_0 = z_0 * sin(t(0));
    double x_1 = z_1 * cos(t(1));
    double y_1 = z_1 * sin(t(1));
    Eigen::Matrix<double, 2, 2> pts;
    of << "plot \"-\" w l lt 5 lw 3" << endl
       << x_0 << " " << y_0 << endl
       << x_1 << " " << y_1 << endl
       << "e" << endl;
}

void Visualizer::save() { of.close(); }

#endif // VISUALIZER_HPP
