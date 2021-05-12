# CellMap

An R pacakge to estimate the cell type proportions of mixture bulk RNA based on pre-computed cell type profiles from sc/sn RNAseq data.

Three main functions are provided within the R package:

  - cellmap: Estimate the cell type proportions of mixture bulk expression based on pre-trained cell type profiles;
  - cellmapTraining: Obtain the cell type profiles of your interests from multiple sc/sn RNAseq datasets by a training process;
  - cellmapOne: Estimate the cell type proportions of mixture bulk based on one reference sc/sn RNAseq dataset.
  
**CellMap is intended to be used for research only and Biogen makes no representation or warranty as to the use or outcome of CellMap**

# Installation
```

# install the cellmap package
devtools::install_github('interactivereport/CellMap')

```

# Pre-build profiles
There are a few pre-build profiles:
  - Major9: Astrocytes, Cardiomyocytes, Endothelial, Hepatocytes, Macrophage, Neuron, Oligodendrocytes, Pancreatic, Skeletal
  - CNS6: Astrocytes, Endothelial, Neuron, Microglia, Oligodendrocytes, Pericyte
  - Neuron3: Inhibitory, Excitatory, Nprogenitor
  - NGN2: iPSC, DAY3, DIV7, DIV7after

Due to the file sizes, all profiles are not installed with the cellmap package. You can obtain those profiles by either directly download from the 'profiles' folder above or call following commands on linux (similar idea for Windows/Mac OS) after installed the cellmap R package:
```

git clone https://github.com/interactivereport/CellMap.git
cd CellMap
./cpProfile.R

```

# manuscript
The manuscript folder contains the modifid version of [***MuSiC***](https://github.com/xuranw/MuSiC), [***SCDC***](https://github.com/meichendong/SCDC) and [***Bisque***](https://github.com/cozygene/bisque), to generate similar output format to the cellmap for comparison/evaluation purpose, as well as parallel implementation of those methods.




