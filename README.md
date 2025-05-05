# InterDILI: drug-induced injury prediction

This model has been trained on a publicly available collection of 5 datasets manually curated for drug-induced-liver-injury (DILI). DILI outcome has been binarised, and ECFP descriptors, together with physicochemical properties have been used to train a random forest classifier which achieves AUROC > 0.9

This model was incorporated on 2024-01-30.

## Information
### Identifiers
- **Ersilia Identifier:** `eos21q7`
- **Slug:** `inter-dili`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Activity prediction`
- **Biomedical Area:** `ADMET`
- **Target Organism:** `Homo sapiens`
- **Tags:** `Toxicity`, `Human`, `Metabolism`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1`
- **Output Consistency:** `Fixed`
- **Interpretation:** Probability of Drug-Induced Liver Injury (DILI), higher score indicates higer risk

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| dili_probability | float | high | Predicted probability of a compound producing drug induced liver injury (DILI) |


### Source and Deployment
- **Source:** `Local`
- **Source Type:** `Replicated`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos21q7](https://hub.docker.com/r/ersiliaos/eos21q7)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos21q7.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos21q7.zip)

### Resource Consumption
- **Model Size (Mb):** `12`
- **Environment Size (Mb):** `783`
- **Image Size (Mb):** `709.72`

**Computational Performance (seconds):**
- 4 inputs: `32.63`
- 20 inputs: `22.28`
- 100 inputs: `22.56`

### References
- **Source Code**: [https://github.com/bmil-jnu/InterDILI](https://github.com/bmil-jnu/InterDILI)
- **Publication**: [https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00796-8](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00796-8)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2024`
- **Ersilia Contributor:** [leilayesufu](https://github.com/leilayesufu)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [None](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos21q7
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos21q7
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
