# reprogramming_ood
Investigation of out-of-distribution problem when aging clocks predict cell reprogramming or embryogenesis.

![Python](https://img.shields.io/badge/python-v3.9+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![made with love](https://img.shields.io/badge/made%20with%20%E2%9D%A4%EF%B8%8F-8A2BE2)


# Title of the paper

A repository containing the code accompanying the research paper "Coming Soon" by Kriukov D., Kuzmina E., ...
(link to biorxiv)

## Description

This repository contains code from a research paper focused on the analysis of neural activity data. Our work delves into the exploration and characterization of *rotational dynamic* prevalent in various neural datasets. We introduce a mathematical framework, informed by our research, designed to assess and quantify the "rotationess" of datasets. This framework leverages **Gyration Numbers**, a complex-valued metric derived from the eigenvalue decomposition of the differential covariance matrix of the data. The resulting **Gyration Plane** facilitates the comparison and simultaneous analysis of multiple datasets.

(image here?)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/ComputationalAgingLab/reprogramming_ood.git
```

2. Navigate to the repository directory:
```bash
cd reprogramming_ood
```

3. (optional) We recommend to reproduce the results using conda environment but you are free to install all necessary libraries manually. In case you choose conda, install it with the following command:

```bash
conda create -n ood python=3.9
```

Then activate the environment:
```bash
conda activate ood
```

4. Install the required packages within the environment:
```bash
pip install -r requirements.txt
```

## Datasets Downloading

Run script that downloads and unzip data: `bash prepare_data.sh`. Or download archive with datasets manually from [here](https://drive.google.com/file/d/11xwxb_m62FymwUeO1vC0KafZ7mog0_wq/view?usp=drive_link). You can also download datasets from GEO (See `paper/Supplementary Table 1.xlsx` for sources of separate datasets).
`dataset.zip` contains preprocessed data, as it is available at origin sources in GEO. 

## Usage
### Data Preparation
Place your data in `./datasets/` directory or execture script `prepare_data.sh` to download data used in the study.
If you would like to work with datasets used in study, you can see how to pre-process and save datasets to `h5` file in notebook `datasets_analysis.ipynb`. If using other dataset, you can use parent class  `NeuralDataset` (from `utils/datasets.py`) and add specific methods to load data and pre-process.

### Content of Repository

- `datasets_analysis.ipynb` - contains tutorial on how to use special class created for working with datasets. It is easy modifiable and allows to use all visualization functions that were used in the study. Have code that was used to render **Fig. S2-S3-S4** and **Fig. 4** from our paper.



### Uncertainty



## Citation
If you use the code or the findings from our paper, **please cite**:

*Plain ref here*

*Bibtex here*

## Contact
For any questions or clarifications, please reach out to: *dmitrii.kriukov@skoltech.ru*

## Datasets Used in Study


## Acknowledgments

Special thanks to Leonid Peshkin for their valuable feedback and suggestions.


## Contributing
We welcome contributions to this repository. If you're found errors in code or experiments, please open an issue to discuss your ideas.


## License
This project is licensed under the MIT License - see the LICENSE file for details.