FROM bentoml/model-server:0.11.0-py38
MAINTAINER ersilia

RUN conda install -c conda-forge openjdk==17.0.11
RUN pip install rdkit==2024.03.3
RUN pip install scikit-learn==1.0.2
RUN pip install numpy==1.21.5
RUN pip install pandas==1.1.3

WORKDIR /repo
COPY . /repo
