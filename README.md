# Predict bioactivity against Main Protease of SARS-CoV-2

MProPred predicts the efficacy of compounds against the main protease of SARS-CoV-2, which is a promising drug target since it processes polyproteins of SARS-CoV-2. This model uses PaDEL-Descriptor to calculate molecular descriptors of compounds. It is based on a dataset of 758 compounds that have inhibition efficacy against the Main Protease, as published in peer-reviewed journals between January, 2020 and August, 2021. Input compounds are compared to compounds in the dataset to measure molecular similarity with MACCS.

This model was incorporated on 2024-07-01.


## Information
### Identifiers
- **Ersilia Identifier:** `eos3nn9`
- **Slug:** `mpro-covid19`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Activity prediction`
- **Biomedical Area:** `COVID-19`
- **Target Organism:** `SARS-CoV-2`
- **Tags:** `COVID19`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1`
- **Output Consistency:** `Fixed`
- **Interpretation:** Gives the pIC50 values for each compound to compare their bioactivity against the main protease

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| mpro_pic50 | float | high | Predicted IC50 value for the inhibition of the Mpro protein of SARS-CoV-2 |


### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos3nn9](https://hub.docker.com/r/ersiliaos/eos3nn9)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos3nn9.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos3nn9.zip)

### Resource Consumption
- **Model Size (Mb):** `83`
- **Environment Size (Mb):** `912`
- **Image Size (Mb):** `1120.65`

**Computational Performance (seconds):**
- 10 inputs: `33.87`
- 100 inputs: `38.25`
- 10000 inputs: `1078.97`

### References
- **Source Code**: [https://github.com/Nadimfrds/Mpropred](https://github.com/Nadimfrds/Mpropred)
- **Publication**: [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10289339/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10289339/)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2023`
- **Ersilia Contributor:** [HarmonySosa](https://github.com/HarmonySosa)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos3nn9
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos3nn9
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
