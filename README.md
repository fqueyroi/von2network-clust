# Source Code and Datasets for submission "Toward Random Walk Based Clustering of Variable-Order Networks" (ICDM 2021)

**Note:** This folder does not contain a main script for all experiments but different ones explained below 

## Dependencies and Setup

1. Python, version >=3 (experiments made with version 3.6.9)
2. [HeapDict](https://pypi.org/project/HeapDict/) python library (command `pip install HeapDict`) (experiments made with version 1.0.1)
3. [Infomap](https://www.mapequation.org/). Users must have a console command `infomap` (experiments made with version 1.3.0). 
4. [LFR Benchmark](https://sites.google.com/site/andrealancichinetti/files). Source code is given in folder 'LFRBenchmark'.

In order to run test experiment reported in Section III, the user must go in folder 'LFRBenchmark' and run command `make` 
See the `ReadMe.txt` in the folder 'LFRBenchmark' for more details.

**The experiments' scripts also use third party codes.**
Files `BuildRulesFast.py` and `BuildNetwork.py` are used for the generation of relevant subsequences in an input dataset and the generation of the corresponding Von networks. The scripts are minor modifications of the ones available at https://github.com/xyjprc/hon (last check June 2021).
LFR Benchmark scripts are used LFR tests cases clustering. We corrected a bug occuring while compiling with `make` (and gcc version 9.3.0) caused by non-void functions in `LFRBenchmark/Sources/print.cpp` and `LFRBenchmark/Sources/histograms.cpp`. Original scripts are available at https://sites.google.com/site/andrealancichinetti/files/binary_networks.tar.gz .

## Datasets

The three datasets used in the paper are available here. An explanation of each dataset is given in Section VI of the paper.

1. `maritime_sequences.csv`: Maritime sequences dataset (default dataset used in the script below)
2. `2011Q1_SEQ.zip` Airport sequences dataset (compressed). 
3. `trajectories_PoliceStation.zip` Taxis sequences dataset (compressed). 

The structure of input files is described at the end of the document. 

## Reproducing the experiments

### LFR Tests cases clustering (Results in Section IV)

In order to reproduce the results given in Fig. 3 (page 5), run

        python3 TestCasesClustering.py

The NMI similarity values used to make Fig. 3 correspond to the following columns in the output:
- 10th Col. (nmi_2o_ns)   : 2-Von input where different codes are assigned to representations of a same location within a given cluster
- 14th Col. (nmi_2o_unif) : 2-Von input where a unique code is assigned to representations of a same location within a given cluster
- 18th Col. (nmi_agg)     : Min 2-Von       input
- 22th Col.,(nmi_fon)     : 2-Fon           input

The boxplots are made using the ggplot2 package of the [R library](https://cran.r-project.org/)

**Warning:** The script creates some temporary files that are not removed at the end.
However, they are written over during each execution.

### Models Accuracy (Results in Section VII)

In order to reproduce the results given in col. *'Acc +- 2sd'* of Table I (page 8), run 

        python3 HONModelsAccuracy.py

To change the dataset used, open file `HONModelsAccuracy.py` and change the variable `filename`.
Using the default value will launch the experiments on the Maritime dataset (file `maritime_sequences.csv`). 
All results are printed inside the Python console.

### Networks and Clustering Comparison (Results in Section VII)

In order to reproduce the results given in Fig. 4, Table I and Table II (page 8), run 

        python3 HONModelsClustering.py

To change the dataset used, open file `HONModelsAccuracy.py` and change the variable `filename`.
Using the default value will launch the experiments on the Maritime dataset (file `maritime_sequences.csv`). 

Network specific outputs are printed inside the Python console.
Node specific statistics (e.g. number of cluster per location) are printed in file `./clusters_stats.csv`.
The reported variables are `id,NET1_nd,NET1_cd,NET2,..` where `NETi` is either the Von2, Agg Von2 or Fon2 network.
- `id`     : id of the location
- `NETi_nd` : number of representations of the location in network `NETi`
- `NETi_nc` : number of clusters found for the location in network `NETi`

The cumulative plots in Fig. 4 are made using the ggplot2 package of the [R library](https://cran.r-project.org/)

The clusters found for each location are printed in file  `./clusters.csv`.
The reported variables are `id;NETi` where `NETi` is either the Von2, Agg Von2 or Fon2 network.
- `id`     : id of the location
- `NETi`   : list of ids of clusters (int) the location belongs to in network `NETi` 

## Using different Datasets of Sequences

The file containing the sequences have the following format

        ID1 L1 L2 L3 ...
        ID2 L5 L6 L2 ...
        ID3 L2 L4 ...
        ...

The first element of each line is the ID of the sequence. 
The rest of the line is the sequence of successive visited locations `Lx`, any string can be used to identified locations.
The separating character (variable `sep`) can be changed in each experiment script.

