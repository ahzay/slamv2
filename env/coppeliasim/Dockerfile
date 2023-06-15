# Use Ubuntu as the base image
FROM ubuntu:latest

# Install CoppeliaSim dependencies
RUN apt-get update -q && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    wget \
    libgl1-mesa-glx \
    libgl1-mesa-dri \
    x11-utils \
    libx11-6 \
    libxcb1 \
    libxau6 \
    libgl1-mesa-dev \
    xvfb \
    dbus-x11 \
    libxkbcommon-x11-0 \
    libavcodec-dev \
    libavformat-dev \
    libswscale-dev \
    libstdc++6 \
    libgcc1 \
    libc6 \
    libicu-dev \
    xz-utils \
    libgl1-mesa-dev \
    && rm -rf /var/lib/apt/lists/*

# Download and install CoppeliaSim
RUN wget https://www.coppeliarobotics.com/files/CoppeliaSim_Edu_V4_5_1_rev4_Ubuntu22_04.tar.xz && \
    tar -xvf CoppeliaSim_Edu_V4_5_1_rev4_Ubuntu22_04.tar.xz && \
    rm CoppeliaSim_Edu_V4_5_1_rev4_Ubuntu22_04.tar.xz

# Set environment variables
ENV COPPELIASIM_ROOT_DIR /CoppeliaSim_Edu_V4_5_1_rev4_Ubuntu22_04
ENV PATH $COPPELIASIM_ROOT_DIR:$PATH
ENV LD_LIBRARY_PATH $COPPELIASIM_ROOT_DIR:$LD_LIBRARY_PATH
ENV QT_DEBUG_PLUGINS=1

