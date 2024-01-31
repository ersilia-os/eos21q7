FROM bentoml/model-server:0.11.0-py38
MAINTAINER ersilia

RUN pip install rdkit==2023.9.2
RUN pip install joblib==1.2.0
RUN pip install numpy==1.24.3
RUN pip install pandas==2.0.3
RUN pip install scikit-learn==1.3.2


WORKDIR /repo
COPY . /repo
