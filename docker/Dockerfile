FROM python:3.9-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN python -m pip install --upgrade pip --no-cache-dir

RUN pip3 install fire==0.5.0 \
	opencv-python-headless==4.6.0.66 \
	scikit-image==0.22.0 \
	dask==2024.2.1 \
	pandas==1.5.3 \
	shapely==2.0.3 \
	imagecodecs==2024.1.1 \
	scanpy==1.9.0 \
	leidenalg==0.10.2 \
	numba==0.59.0 \
	'matplotlib<3.7'
