# Hamming Weight Tree (HWT)
Implementation of the Hamming Weight Tree in C++.

HWT supports two fundumental queries in binary datasets: 1) K nearest neighbor search 2) Insertion of new binary codes

### Compilation
Compiling and building the project requires cmake, make, hdf5 and hdf5-dev. To compile go to hwt folder and run cmake:

```
cmake CMakeLists.txt
make
```
This should create two executable files 1) *hwt*, for running the Hamming weight tree, and 2) *linscan* for running linear scan

### Usage
##### hwt
hwt provides an implementaiton of the Hamming weight tree which supports fast Hamming nearest neighbor search as well as fast insertion of new items to the data strucutre
