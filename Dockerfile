# FROM nvcr.io/nvidia/pytorch:20.10-py3
FROM continuumio/miniconda3

WORKDIR /app

COPY environment.yml .

RUN conda env create -f environment.yml

SHELL ["conda", "run", "-n", "collectnet", "/bin/bash", "-c"]
# RUN conda env create -f environment.yml --verbose
# RUN conda clean -afy

# SHELL ["conda", "run", "-n", "collectnet", "/bin/bash", "-c"]

COPY . .

CMD ["Rscript", "infer.r", "default_value1", "default_value2", "default_value3"]