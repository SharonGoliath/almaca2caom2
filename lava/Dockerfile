FROM opencadc/astropy:3.6-slim

WORKDIR /home/da_data/lava-RHEL7/LAVA-RHEL7

COPY requirements.txt .

RUN pip install -r requirements.txt

RUN apt-get update -y && apt-get install -y vim

RUN pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools && \
    pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatasks
