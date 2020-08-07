FROM ubuntu:bionic
RUN apt update && apt install -y wget vim gnupg2 software-properties-common
RUN wget -O - https://dl.openfoam.org/gpg.key | apt-key add --
RUN add-apt-repository http://dl.openfoam.org/ubuntu
RUN apt install -y openfoam7 && apt clean
