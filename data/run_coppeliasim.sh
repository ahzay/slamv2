#/bin/bash 
xhost +local:root
mkdir -p $(pwd)/shared_data/scans
docker run --ulimit nofile=1024:524288 -it --rm   --env DISPLAY=$DISPLAY  --privileged  --volume /tmp/.X11-unix:/tmp/.X11-unix   --volume /home/user/.Xauthority:/root/.Xauthority   --network=host   --env ROS_MASTER_URI=http://localhost:11311   --env ROS_IP=localhost   --env ROS_HOSTNAME=localhost  -v "/$(pwd)/shared_data":"/data"  -e RMW_IMPLEMENTATION=rmw_fastrtps_cpp -e ROS_DOMAIN_ID=5 -e QT_QPA_PLATFORM=xcb coppeliasim coppeliaSim.sh $(ls -t /data/*ttt | head -n1)
mkdir -p ../src/cmake-build-debug-docker/scans
cp shared_data/scans/* ../src/cmake-build-debug-docker/scans
