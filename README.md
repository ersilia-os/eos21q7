# InterDILI: drug-induced injury prediction

Flags the risk that a compound causes drug-induced liver injury, the leading cause of safety-driven withdrawals from the market. Lee and Yoo built InterDILI from substructure and physicochemical descriptors with an attention mechanism, reporting AUROC between 0.88 and 0.97 across five public DILI datasets. Ersilia replicates it with a LazyQSAR classifier trained on the merged set of 1,850 compounds, 952 of them hepatotoxic, with a five-fold out-of-fold AUROC of 0.86.

This model was incorporated on 2024-01-30.
Last packaged on 2026-11-23.

## Information
### Identifiers
- **Ersilia Identifier:** `eos21q7`
- **Slug:** `inter-dili`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Activity prediction`
- **Biomedical Area:** `ADMET`
- **Target Organism:** `Homo sapiens`
- **Tags:** `Toxicity`, `Metabolism`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1`
- **Output Consistency:** `Fixed`
- **Interpretation:** Rank score calibrated against the Ersilia reference library of 50K compounds, where 0.65 beats 99% of it. Higher values mean a greater risk that the compound causes drug-induced liver injury.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| dili_rank | float | high | Rank score of how likely the compound is to cause drug-induced liver injury (DILI). Calibrated against 50K compounds from the Ersilia reference library so that a score of 0.65 is better than 99% of them. Recommended threshold: 0.65. |


### Source and Deployment
- **Source:** `Local`
- **Source Type:** `Replicated`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos21q7](https://hub.docker.com/r/ersiliaos/eos21q7)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos21q7.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos21q7.zip)
- **Apptainer**: [https://models-sif.s3.eu-north-1.amazonaws.com/eos21q7_v3.sif](https://models-sif.s3.eu-north-1.amazonaws.com/eos21q7_v3.sif)

### Resource Consumption
- **Model Size (Mb):** `31`
- **Environment Size (Mb):** `6954`
- **Image Size (Mb):** `3288.63`

**Computational Performance (seconds):**
- 10 inputs: `35.8`
- 100 inputs: `27.4`
- 10000 inputs: `181.22`

### References
- **Source Code**: [https://github.com/bmil-jnu/InterDILI](https://github.com/bmil-jnu/InterDILI)
- **Publication**: [https://doi.org/10.1186/s13321-023-00796-8](https://doi.org/10.1186/s13321-023-00796-8)
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
