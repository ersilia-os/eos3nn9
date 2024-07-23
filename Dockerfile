FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN conda install -c conda-forge openjdk
RUN pip install rdkit
RUN pip install pandas==1.1.3
RUN pip install numpy==1.21.5
RUN pip install scikit-learn==1.0.2
RUN pip install isaura==0.1
RUN pip install install-jdk

ENV PATH="$PATH:bin:/root/bin:jvm/bin:/root/jvm/bin:miniconda3/envs/ersilia/lib/jvm/bin/java:/root/miniconda3/envs/ersilia/lib/jvm/bin/java:miniconda3/pkgs/openjdk-22.0.1-hb622114_0/lib/jvm/bin:/root/miniconda3/pkgs/openjdk-22.0.1-hb622114_0/lib/jvm/bin"

WORKDIR /repo
COPY . /repo
