# Cell Based Comparison 2017

This repository contains the code used to produce results in the following paper:
[Comparing individual-based approaches to modelling the self-organization of multicellular tissues](https://doi.org/10.1371/journal.pcbi.1005387)

The authors on this paper are:
 - James M. Osborne
 - Alexander G. Fletcher
 - Joe M. Pitt-Francis
 - Philip K. Maini
 - David J. Gavaghan


## Code to reproduce simulations

The following test files run simulations from the paper:

 - [TestCellSortingLiteratePaper.hpp](./test/TestCellSortingLiteratePaper.hpp)
 - [TestCylindricalCryptLiteratePaper.hpp](./test/TestCylindricalCryptLiteratePaper.hpp)
 - [TestDeltaNotchLiteratePaper.hpp](./test/TestDeltaNotchLiteratePaper.hpp)
 - [TestMorphogenMonolayerLiteratePaper.hpp](./test/TestMorphogenMonolayerLiteratePaper.hpp)


## Continuous testing

This project is regularly compiled, and short versions of each simulation are regularly run every time the main Chaste repository is updated.
If you encounter a problem running the code in this project, please [open an issue](https://github.com/Chaste/CellBasedComparison2017/issues).

The original parameter values are available in a block near the top of each file.


## Running these simulations

You must first compile [Chaste](https://github.com/Chaste/Chaste).
Visit our [getting started](https://chaste.github.io/getting-started/) page for full details.

From the `projects` directory, clone this repository:

```bash
projects$ git clone https://github.com/Chaste/CellBasedComparison2017.git
```

From your Chaste build directory, configure, compile and run with any configuration you like, for instance:

```bash
build$ cmake -DCMAKE_BUILD_TYPE=Release -DChaste_ERROR_ON_WARNING=OFF /path/to/chaste
build$ make -j4 project_CellBasedComparison2017
build$ ctest -j4 -L project_CellBasedComparison2017
```
