FROM ubuntu:24.10

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
COPY ./model/ /home/copied_stuff

#install Intel stuff
RUN mkdir -p $CRYOSHOME
WORKDIR $CRYOSHOME
RUN cp -r /home/copied_stuff/intel . && sh ./intel/intel-oneapi-hpc-toolkit-2025.0.1.47_offline.sh -a --silent --cli --eula accept && sh ./intel/intel-oneapi-base-toolkit-2025.0.1.46_offline.sh -a --silent --cli --eula accept && rm -r intel/
#&& . /opt/intel/oneapi/setvars.sh
RUN cp -r /home/copied_stuff/les-lsm-large_domain/* $CRYOSHOME

#compile the code
RUN apt install -y makedepf90


#create directory where palm_exchange folder gets mapped
#RUN mkdir $INSTALLFOLDER/JOBS
#WORKDIR $INSTALLFOLDER
#RUN mkdir -p $INSTALLFOLDER/JOBS/example_cbl && cp -r /home/palmuser/jobtocopy/* $INSTALLFOLDER/JOBS/example_cbl
#cp -r $PALMFOLDER/packages/palm/model/tests/cases/example_cbl $INSTALLFOLDER/JOBS

#do PALM stuff
#RUN $INSTALLFOLDER/bin/palmrun -b -a "d3#" -r example_cbl -q medium -X 256 -T 64 -t 1800 -m 60000 -v -F
#RUN $INSTALLFOLDER/bin/palmrun -r example_cbl -c default -m 60000 -t 1800 -q medium -i 16683 -U root -a "d3#" -G "PALM 24.04" -X 2 -T 2 -B -v
#RUN mv $INSTALLFOLDER/tmp/example_cbl* $INSTALLFOLDER/tmp/job_cscs
#WORKDIR $INSTALLFOLDER/tmp/job_cscs
#RUN rm -r *.f90 && rm -r *.mod && rm -r *.o && rm -r *.c && rm -r *.tar
#RUN chmod +x jobfile.*
