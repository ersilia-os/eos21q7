# interpretable prediction of drug-induced liver injury through permutation feature importance and attention mechanism

This study focused on interpreting the model prediction to analyze which features could potentially cause DILI and indicated that the model not only captured the previously known properties to be related to DILI but also proposed a new DILI potential substructural of physicochemical properties.

## Identifiers

* EOS model ID: `eos21q7`
* Slug: `inter_dili`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Classification`
* Output: `Probability`
* Output Type: `Float`
* Output Shape: `Single`
* Interpretation: The model predicts the probability of Drug-Induced Liver Injury (DILI) based on compound features, providing a probability score ranging from 0 to 1, with higher values indicating greater risk

## References

* [Publication](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00796-8)
* [Source Code](https://github.com/bmil-jnu/InterDILI/tree/main)
* Ersilia contributor: [leilayesufu](https://github.com/leilayesufu)

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos21q7)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos21q7.zip)
* [DockerHub](https://hub.docker.com/r/ersiliaos/eos21q7) (AMD64, ARM64)

## Citation

If you use this model, please cite the [original authors](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00796-8) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a None license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!