FROM python:3-slim
MAINTAINER ADAM GREEN

COPY requirements.txt requirements.txt

RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    libhdf5-dev \
    libblas-dev \
    liblapack-dev \
    libxft-dev \
    wget \
    vim \
    libsm6 \
    libxext6 \
    libxrender-dev &&\
    wget https://github.com/opensourcedesign/fonts/raw/master/gnu-freefont_freesans/FreeSansBold.ttf &&\
    mv FreeSansBold.ttf /usr/share/fonts/truetype &&\
    rm -rf /var/lib/apt/lists/* &&\
    pip3 install --upgrade pip &&\
    pip3 install -r requirements.txt


COPY init.sh init.sh

WORKDIR /work
#RUN export PATH=$PATH:/app/src/data && echo "export PATH=$PATH" >> /etc/bash.bashrc


CMD ["bash","call_defect.sh"]
#CMD ./call_defect
#CMD bash
