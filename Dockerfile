FROM rocker/r-ver:4.1.0

RUN apt-get update && apt-get install -y wget bzip2 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    /opt/conda/bin/conda init

ENV PATH /opt/conda/bin:$PATH

COPY environment.yml /app/environment.yml

RUN conda env create -f /app/environment.yml

RUN echo "conda activate scrna" >> ~/.bashrc
ENV PATH /opt/conda/envs/scrna/bin:$PATH

COPY . /app

WORKDIR /app

CMD ["Rscript", "inference.r"]