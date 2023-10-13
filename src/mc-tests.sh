#/bin/bash
./cmake-build-debug-docker/mc-simulator -5 5 -5 5 240 684 ./simulated > /dev/null
./cmake-build-debug-docker/mappernode 1 0 50 1 0.25 0.25 1 1 100 5 3 ./simulated > /dev/null
./cmake-build-debug-docker/evaluator
