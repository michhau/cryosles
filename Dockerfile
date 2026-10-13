FROM ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive
ENV CRYOSHOME="/home/cryosles"
#ENV PALMFOLDER="/home/palmuser/palm/current_version/palm_model_system-v24.04/"
#ENV INSTALLFOLDER="/home/palmuser/palm/current_version/"
#ENV FFLAGS="-fallow-argument-mismatch"

RUN apt-get -y update && apt-get install -y apt-utils

#install dependencies
RUN apt-get install -y cmake pkg-config build-essential wget vim
#apt-get install -y gfortran g++ make cmake coreutils libnetcdff-dev netcdf-bin libfftw3-dev python3-pip python3-pyqt5 flex bison build-essential python3.12-venv python3-numpy wget
#RUN apt-get install -y openmpi-bin libopenmpi-dev

#copy stuff
WORKDIR /home
RUN mkdir -p /home/copied_stuff
COPY ./libs/ /home/copied_stuff/
COPY ./model /home/copied_stuff

#install Intel stuff
RUN mkdir -p $CRYOSHOME
WORKDIR /home/copied_stuff
#RUN sh ./intel/intel-oneapi-hpc-toolkit-2026.1.1.33_offline.sh -a --silent --cli --eula accept && 
RUN sh ./intel/intel-oneapi-toolkit-2026.1.1.33_offline.sh -a --silent --cli --eula accept && rm -r /home/copied_stuff/intel/
#RUN . /opt/intel/oneapi/setvars.sh

#build model
WORKDIR $CRYOSHOME
RUN cp -r /home/copied_stuff/* $CRYOSHOME

#compile the code
RUN apt install -y makedepf90
