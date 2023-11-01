//
// Created by user on 10/12/23.
//
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <algorithm>
#include "../mappernode/tools.h"

using namespace std;
using namespace Eigen;

struct Point {
    double x;
    double y;
};

struct Circle {
    Point center;
    double radius;
};

struct Rectangle {
    Point bottom_left;
    Point top_right;
};

Rectangle findSmallestRectangle(const Circle &c1, const Circle &c2) {
    Rectangle rect;
    rect.bottom_left.x = min(c1.center.x - c1.radius, c2.center.x - c2.radius);
    rect.bottom_left.y = min(c1.center.y - c1.radius, c2.center.y - c2.radius);
    rect.top_right.x = max(c1.center.x + c1.radius, c2.center.x + c2.radius);
    rect.top_right.y = max(c1.center.y + c1.radius, c2.center.y + c2.radius);
    return rect;
}


Vector<double, 6> load_params(string fn) {
    ifstream f(fn);
    if (!f.is_open())
        throw runtime_error("failed to open file");
    Vector<double, 6> v;
    for (int i = 0; i < 6; i++)
        f >> v(i);
    return v;
}

Circle se_to_circle(Vector<double, 6> p) {
    return {
            Point({p(0), p(1)}),
            max(p(3), p(4))
    };
}

bool is_inside_se(const Vector<double, 6> &p, Point t) {
    double x = t.x;
    double y = t.y;
    double f1 = ((x - p(0)) * cos(p(2)) + (y - p(1)) * sin(p(2))) / p(3);
    double f2 = ((x - p(0)) * sin(p(2)) - (y - p(1)) * cos(p(2))) / p(4);
    double buf = pow(pow(f1, 2.), 1. / p(5)) + pow(pow(f2, 2.), 1. / p(5));
    double e = sqrt(p(3) * p(4)) * (pow(buf, p(5)) - 1.);
    if (e <= 0)
        return true;
    else return false;
}

double evaluate() {
    Vector<double, 6> est, tru;
    tru = load_params("ground_truth.txt");
    est = load_params("ground_estimate.txt");
    Circle c_est = se_to_circle(est),
            c_tru = se_to_circle(tru);
    Rectangle r = findSmallestRectangle(c_est, c_tru);
    double precision = 1e-2;
    int inside_tru = 0;
    int inside_both = 0;
    for (double x = r.bottom_left.x; x <= r.top_right.x; x += precision) {
        for (double y = r.bottom_left.y; y <= r.top_right.y; y += precision) {
            Point p({x, y});
            if (is_inside_se(tru, p)) {
                inside_tru++;
                if (is_inside_se(est, p))
                    inside_both++;
            }
        }
    }
    return fmod(((double) inside_both / (double) inside_tru), 1.);
}

int main(int argc, char *argv[]) {
    int num = 1000;
    //double avg_accuracy = 0;
    double sum = 0;
    double sum_of_squares = 0;
    int valid_entries = 0;
    for (int i = 0; i < num; i++) {
        system("./mc-simulator -5 5 -5 5 240 684 ./simulated > /dev/null");
        system("./mappernode 1 0 50 1 100 100 100 1 100 5 1 $(nproc) ./simulated > /dev/null");
        system("mv 0.png simulated/$(date +%s).png");
        double accuracy = evaluate();
        cout << "accuracy is: " << accuracy << endl;
        if (accuracy != 0 && isfinite(accuracy)) {
            sum += accuracy;
            sum_of_squares += pow(accuracy, 2);
            valid_entries++;

            double avg_accuracy = sum / valid_entries;
            double variance = (sum_of_squares - pow(sum, 2) / valid_entries) / valid_entries;
            double std_dev = sqrt(variance);
            double std_error = std_dev / sqrt(valid_entries);

            cout << "average accuracy is: " << avg_accuracy << endl;
            cout << "standard error is: " << std_error << endl;
        }
    }
}