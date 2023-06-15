FROM docker.io/archlinux/archlinux:latest
ENV MAKEFLAGS="-j$(nproc)"
#RUN pacman -Syu --noconfirm 
# Install base-devel and git for AUR helper
RUN pacman -Syu --noconfirm base-devel git parallel gnuplot gdb

# Create a new user for the AUR helper
RUN useradd -m builder && \
    echo "builder ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/builder

USER builder
WORKDIR /home/builder

# Install an AUR helper
RUN git clone https://aur.archlinux.org/yay.git && \
    cd yay && \
    makepkg -si --noconfirm

# Install the required AUR packages
RUN yay -S --noconfirm autodiff

USER root

# Install the required official repository packages
RUN pacman -S --noconfirm cmake boost ceres-solver python python-numpy gcc openmp
RUN pacman -S --noconfirm qtcreator 

CMD ["bash"]
