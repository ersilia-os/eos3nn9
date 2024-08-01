# Predict bioactivity against Main Protease of SARS-CoV-2

MProPred predicts the efficacy of compounds against the main protease of SARS-CoV-2, which is a promising drug target since it processes polyproteins of SARS-CoV-2. This model uses PaDEL-Descriptor to calculate molecular descriptors of compounds. It is based on a dataset of 758 compounds that have inhibition efficacy against the Main Protease, as published in peer-reviewed journals between January, 2020 and August, 2021. Input compounds are compared to compounds in the dataset to measure molecular similarity with MACCS.

## Identifiers

* EOS model ID: `eos3nn9`
* Slug: `mpro-covid19`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Regression`
* Output: `Score`
* Output Type: `Float`
* Output Shape: `List`
* Interpretation: Gives the pIC50 values for each compound to compare their bioactivity against the main protease

## References

* [Publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10289339/)
* [Source Code](https://github.com/Nadimfrds/Mpropred)
* Ersilia contributor: [HarmonySosa](https://github.com/HarmonySosa)

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos3nn9)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos3nn9.zip)
* [DockerHub](https://hub.docker.com/r/ersiliaos/eos3nn9) (AMD64, ARM64)

## Citation

If you use this model, please cite the [original authors](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10289339/) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a MIT license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!