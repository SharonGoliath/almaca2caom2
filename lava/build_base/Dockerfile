FROM python:3.6-slim-jessie

WORKDIR /home/da_data/lava-RHEL7/LAVA-RHEL7

ENV ENV="PROD"

COPY entrypoint.sh /entrypoint.sh

RUN apt-get update \
    && apt-get install -y gcc libc-dev libffi-dev libgfortran3 libmagic-dev libssl-dev libxslt-dev make vim \
    && apt-get clean \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip to the latest.
RUN pip install --upgrade pip
RUN pip install --no-cache-dir "astropy>=3.2.3,<4.0" \
    && pip install --no-cache-dir "regions"

  
COPY requirements.txt .

RUN pip install -r requirements.txt

RUN pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools && \
    pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatasks

RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
