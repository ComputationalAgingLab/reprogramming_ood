# reprogramming_ood
Investigation of out-of-distribution problem when aging clocks predict cell reprogramming or embryogenesis.

![Python](https://img.shields.io/badge/python-v3.9+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![made with love](https://img.shields.io/badge/made%20with%20%E2%9D%A4%EF%B8%8F-8A2BE2)


# Title of the paper

A repository containing the code accompanying the research paper "On the Rotational Structure in Neural Data" by Kuzmina E., Kriukov D., Lebedev M.
[bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.11.557230v2)

## Description

This repository contains code from a research paper focused on the analysis of neural activity data. Our work delves into the exploration and characterization of *rotational dynamic* prevalent in various neural datasets. We introduce a mathematical framework, informed by our research, designed to assess and quantify the "rotationess" of datasets. This framework leverages **Gyration Numbers**, a complex-valued metric derived from the eigenvalue decomposition of the differential covariance matrix of the data. The resulting **Gyration Plane** facilitates the comparison and simultaneous analysis of multiple datasets.

<p align="center">
<img src="image/intro.png" alt>

</p>
<p align="center">
<em>The brain areas that were explicitly studied with rotational dynamics approach in rhesus monkeys, humans and rodents. </em>
</p>


## Installation

1. Clone the repository:
```bash
git clone https://github.com/ComputationalAgingLab/reprogramming_ood.git
```

2. Navigate to the repository directory:
```bash
cd reprogramming_ood
```

3. (optional) We recommend ...

```bash
conda create -n ood python=3.9
```

Then activate ...
```bash
conda activate ood
```


4. Install the required packages:
```bash
pip install -r requirements.txt
```

## Datasets Downloading

1. Run script that downloads and unzip data: `bash prepare_data.sh` or `bash prepare_preproc_datasets.sh`. Or download archive with datasets manually from [here](https://drive.google.com/drive/folders/1AWO8XZpLBW1fkp5ylF6-w6J8gYcnnOkp?usp=sharing). You can also download datasets from the original source from here ... # TODO link to table.
`dataset.zip` contains raw data, as it is available at origin sourse. `preprocessed2h5.zip` contains already prerocessed datasets, saved to .h5 file as dictionary.
 
2. Add path to downloaded datasets to `datasets_config.py` file for convinience.

## Usage
### Data Preparation
Place your data in `./datasets/` directory or execture script `prepare_data.sh` to download data used in the study.
If you would like to work with datasets used in study, you can see how to pre-process and save datasets to `h5` file in notebook `datasets_analysis.ipynb`. If using other dataset, you can use parent class  `NeuralDataset` (from `utils/datasets.py`) and add specific methods to load data and pre-process.

### Content of Repository

- `datasets_analysis.ipynb` - contains tutorial on how to use special class created for working with datasets. It is easy modifiable and allows to use all visualization functions that were used in the study. Have code that was used to render **Fig. S2-S3-S4** and **Fig. 4** from our paper.

- `gyration_plane_tutorial.ipynb` - contains tutorial on how to use and interpret **Gyration Number** concept. How to plot datasets in **Gyration Plan** and compare them. Also, contains bonus explanation of **Curvature** concept.

- `traveling_wave_model.ipnyb` - Have code that was used to render **Fig. 3, 5** and **Fig. S5-S6-S7** from our paper.


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