FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN pip install rdkit
RUN pip install pandas==1.1.3
RUN pip install numpy==1.21.5
RUN pip install scikit-learn==1.0.2
# RUN pip install install-jdk
RUN conda install -c conda-forge openjdk
# ENV PATH="$PATH:bin:./bin:jvm/bin"

WORKDIR /repo
COPY . /repo
