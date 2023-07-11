//
// Created by user on 7/11/23.
//

#include "cmdlineoptions.h"
#include <iostream>

CmdLineOptions::CmdLineOptions(int argc, char **argv) {
    if (argc < 10)
        throw std::runtime_error("Not enough input arguments.");
    scan_number = stoi(argv[1]);
    prune = stoi(argv[2]);
    init_npoints = stoi(argv[3]);
    stride = stoi(argv[4]);
    distance_tolerance = strtof(argv[5], nullptr);
    angle_tolerance = strtof(argv[6], nullptr);
    env_multiplier = strtof(argv[7], nullptr);
    max_scan_distance = strtof(argv[8], nullptr);
    data_folder = string(argv[9]);
}
