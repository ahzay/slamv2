# Geometric SLAM
A novel method for solving the SLAM problem in 2D environments using measurements from a Laser Range Finder (LRF)

## Usage
Docker containers are provided for both the simulation and the build environment. To build all the required containers, make sure Docker is installed and the user is added to the Docker group (with appropriate permissions) then run:
```
cd ./env
./build_all.sh
```

## Notes
- `Invalid MIT-MAGIC-COOKIE-1 keyqt.qpa.xcb: could not connect to display :0` can be solved by running `xhost +local:`


# TODOs
- [ ] complete README.md
