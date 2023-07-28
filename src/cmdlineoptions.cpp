//
// Created by user on 7/11/23.
//

#include "cmdlineoptions.h"
#include <iostream>

CmdLineOptions::CmdLineOptions(int argc, char **argv) {
    if (argc < 14)
        throw std::runtime_error("Not enough input arguments.");
    init_npoints_multiplier = strtof(argv[8], nullptr);
    scan_number = stoi(argv[1]);
    prune = stoi(argv[2]);
    init_npoints = stoi(argv[3]);
    stride = stoi(argv[4]);
    distance_tolerance = strtof(argv[5], nullptr);
    angle_tolerance = strtof(argv[6], nullptr);
    eucledian_tolerance = strtof(argv[7], nullptr);
    max_scan_distance = strtof(argv[9], nullptr);
    data_longevity = stoi(argv[10]);
    npoints_mult = stoi(argv[11]);
    nthreads = stoi(argv[12]);
    data_folder = string(argv[13]);
}
