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
##### Hamming Weight Tree
`hwt` provides an implementaiton of the Hamming weight tree which supports fast Hamming nearest neighbor search as well as fast insertion of new items to the data strucutre

```
./hwt <.mat inout file> <h5 output file> -K <#neighbors> -B <length of code> -Q <#queries> -T <capacity>
```

```
e.g.:/hwt binary.mat out.h5 -K 10 -B 64 -Q 100 -T 1000
```
The input file should be a .mat file containing the matrix of binary codes `B` and the matrix of query points `Q`. In both matrices the columns reporsent the data points and rows represent the bits. To reduce the storage cost, each 8 bits of binary code must be compact into a byte, ranging from 0 to 255. Therefore, the number of rows is equal to (length of code)/8. `K` is the number of nearest neighbors to retreive, `B` is the legnth of codes, `Q` is the number of queries and `T` is the maximum number of nodes that can assigned to a leaf node.

##### Linear Scan
`linscan` provides an implementation of the linear search to find the K nearest neighbors for dynamic datasets.
```
./linscan  <.mat inout file> <h5 output file> -K <#neighbors> -B <length of code> -Q <#queries>

e.g.:/hwt binary.mat out.h5 -K 10 -B 64 -Q 100
```

The syntax is similar to `hwt` except that it does not require to specify the capacity.
