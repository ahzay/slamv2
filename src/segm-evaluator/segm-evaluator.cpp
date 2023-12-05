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

struct Segment {
    vector<Vector2d> ms; // vector of measurements
};

vector<Segment> load_segmentation(string fn) { // index 0 is not associated
    std::ifstream ifs(fn);
    vector<Segment> segments;
    string line;
    int n;
    int z = 0;
    if (ifs.is_open()) {
        ifs >> n;
        segments.resize(z + 1);
        for (int i = 0; i < n; n++) {
            Vector2d m;
            int idx;
            ifs >> m[0] >> m[1] >> idx;
            if (idx != z) {
                z++;
                segments.resize(z + 1);
            }
            segments[z].ms.push_back(m);
        }
    } else throw std::runtime_error("Could not open file: " + fn);
    return segments;
}

int main(int argc, char *argv[]) {
    auto est = load_segmentation("ground_segm_estimate.txt");
    auto tru = load_segmentation("ground_segm_truth.txt");

}