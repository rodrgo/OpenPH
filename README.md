# OpenPH: Persistent Homology computation on GPUs

OpenPH is a library for fast persistent homology computation on GPUs. OpenPH provides a `CUDA-C` implementation of `pms`, which is a super-fast parallel algorithm for reduction of boundary matrices designed for GPUs. Additionally, it includes vanilla implementations of other classical reduction algorithms like `standard`, `twist` and `ph-row`. Please, see the paper on [Arxiv](https://arxiv.org/abs/1708.04710). 

## Installation

### OpenPH

The following programs are required to build OpenPH: `make`, `python`, `nvcc`, `wget`, `matlab`. Unfortunately, `matlab` is currently the only interface for OpenPH (this is due to inherited practices in numerical computing groups). However, future versions will include `python` and `julia` interfaces. Maybe you can help writing them :)

To setup the project and build OpenPH, please run:

```bash
sh install.sh
```

The script will do its best to automatically populate `src/cuda/config` with your system parameters (like location of `CUDA` and `MATLAB` installation, for example). If this fails, please look at `src/config` and fill in the blanks. 

### Numerics

The following programs are required to build the dependencies and run the numerics: `cmake`, `ant`. 

To install the dependencies and run the numerics, please do

```bash
cd numerics
sh install.sh
sh run.sh
```

The following libraries are needed to run the numerical experiments:

* [Javaplex](https://github.com/appliedtopology/javaplex) - To generate pointclouds, simplicial complexes and boundary matrices.
* [Phat](https://bitbucket.org/phat-code/phat) - To provide a baseline.

Additionally, some experiments use the datasets in [PH-roadmap](https://github.com/n-otter/PH-roadmap).

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Citing

If you would like to cite this software, please cite this [Arxiv paper](https://arxiv.org/abs/1708.04710). 

## License

[MIT License](https://github.com/rodrgo/openph/blob/master/LICENSE)

