In this folder there are 6 simulated datasets for user to test the cellmap.
One set of pure (single cell type) and one set of mixture (mix cell types) pseudo bulk samples were created from datasets used to generat profiles (Major9, CNS6 and Neuron3) respectively.
The R object (rds) file contains an pseudo bulk expression matrix, each column name indicates the dataset it was created from. The number of cells used to created each sample is in the corresponding txt file (*.rate, table separated by '\t') with each row is a cell type, each column is a pseudo sample which is consistent with the corresonding pseudo expression matrix.

Sample number:
      Major9      CNS6      Neuron3
Pure   220         220        110
Mix    80           45        45
