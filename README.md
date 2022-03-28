# `haddock2mmcif`

🚧 Attention 🚧

This code is under development and subject to change. Please report any bugs by opening an issue in the repository.

* * *

Encode information from a HADDOCK run into a `.cif` to be deposited in PDB-Dev.

## Currently the follow information is encoded in the `.cif`

- Models are represented as:

    1. Whole structure as rigid
    2. Interface as flexible, defined by the `flcut` parmameter in `run.cns`

- Top4 models from each clusters within `ModelGroup`, ranked by their top4 HADDOCK score
- Restraints as `Dataset` represented as `DerivedDistanceRestraint`

    1. Ambiguous, Active/Passive as `ResidueFeature` with probability defined by `ncvpart` in `run.cns` read from `ambig.tbl`
    2. Unambiguous, the `ResidueFeature` are read from `unambig.tbl`


## To-be-implemented

- HADDOCK-score

* * *

## Installation

```
$ git clone https://github.com/haddocking/haddock2mmcif
$ cd haddock2mmcif
$ python setup.py develop
$ cd tools
$ g++ -o contact-chainID contact-chainid.cpp
$ cd ..
```

## Usage

The input for `haddock2mmcif` is the folder of a run-directory, either downloaded from the web server or executed locally.

```
$ cd example_data
$ tar zxvf 6269-E2A-HPR.tgz
$ haddock2mmcif --output example.cif E2A-HPR/
```
