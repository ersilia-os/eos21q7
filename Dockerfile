FROM bentoml/model-server:0.11.0-py311
MAINTAINER ersilia

RUN pip install rdkit==2023.9.4
RUN pip install joblib==1.3.2
RUN pip install numpy==1.26.4
RUN pip install pandas==2.2.0
RUN pip install scikit-learn==1.4.0


WORKDIR /repo
COPY . /repo
