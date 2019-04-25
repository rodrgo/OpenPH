# OpenPH: Persistent Homology computation on GPUs

OpenPH is a library for fast persistent homology computation on GPUs. OpenPH provides a `CUDA-C` implementation of `pms`, which is a provably convergent parallel algorithm for boundary matrix reduction tailored for GPUs, see the [Arxiv paper](https://arxiv.org/abs/1708.04710). Additionally, it includes vanilla implementations of other classical reduction algorithms like `standard`, `twist` and `ph-row`.

## Installation

### OpenPH

The following programs are required to build OpenPH: `make`, `python`, `nvcc`, `wget`, `matlab`. At this moment, `matlab` is the only interface for OpenPH, but we're working hard to make this project `matlab` independent. In particular, we plan to add `python` and `julia` APIs by release `1.0.0`.

To setup the project and build OpenPH for `matlab`, please run:

```bash
sh build.sh
```

The script will do its best to automatically generate a `config` file (in the `src/matlab` folder) with your system and GPU parameters. If this fails, please look at `src/matlab` and fill in the blanks. 

### Numerics

The following programs are required to build the dependencies and run the numerics: `cmake`, `make`, `ant`, `java`, `openmp`, `epstopdf`.

To install the dependencies and run the numerics, please do

```bash
cd numerics
bash install.sh
bash run.sh
```

The following libraries will be installed as the are needed to run the numerical experiments. The dependencies are listed below each library.

* [Javaplex](https://github.com/appliedtopology/javaplex) - To generate pointclouds, simplicial complexes and boundary matrices.
    * `java`, `ant`
* [Phat](https://bitbucket.org/phat-code/phat) - Benchmarking baseline.
    * `cmake`, `make`
* [DIPHA](https://github.com/DIPHA/dipha) - Benchmarking datasets.
    * `cmake`, `make`, `openmp`

## Overview

The current (vanilla) version of OpenPH uses a sparse representation of the boundary matrix that allocates a fixed number of memory slots per column. This induces a segmentation of the data that can be naturally split across $m$ processors (one column per GPU thread) which permits a fast routine for column addition. Parallelisation is achieved by identifying a set of `pivots` (columns in the matrix that are already reduced) at the start of the algorithm and using them to apply a left-to-right column operation wherever possible. After each left-to-right column operation the algorithm checks whether the column can be added to the set of `pivots`. 

It is shown in the numerical experiments that this strategy yields a very fast algorithm. However, our vanilla version suffers from one main drawback:

1. The number of memory slots per column has to be supplied at the beginning. In release `0.0.1` a value of `p*MAX_NNZ` was chosen, where `MAX_NNZ` is the number of non-zeros of the *least sparse* column in the matrix at the start of the algorithm and `p` is an integer chosen heuristically. 

This is not only sub-optimal from the storage complexity point of view, but can also make the algorithm unstable if the number of non-zeros in a given column grows beyond `p*MAX_NNZ` during reduction. The main objective prior to release `1.0.0` is to improve storge complexity and re-allocate memory adaptively.

## Contributing

Please contribute to the project! :) 

Detailed contributing guidelines will be established in release `1.0.0`. If in the meantime, you wish to contribute, please feel free to do so by following the standard [fork & pull request workflow](https://gist.github.com/Chaser324/ce0505fbed06b947d962). If you intend to submit a contribution, it might be a good idea to open a [draft pull-request](https://github.blog/2019-02-14-introducing-draft-pull-requests/) and [write us an email](mailto:r.mendozasmith@gmail.com) to discuss.

## Citing

```bibtex
@article{mendoza2017parallel,
  title={Parallel multi-scale reduction of persistent homology filtrations},
  author={Mendoza-Smith, Rodrigo and Tanner, Jared},
  journal={arXiv preprint arXiv:1708.04710},
  year={2017}
}
```

## License

We want everyone to use this code and to contribute to this project, but we also want to keep patent trolls away. This is why we picked an Apache license rather than an MIT license. If you have a strong point of view about licensing of numerical computing libraries, please [write us an email](mailto:r.mendozasmith@gmail.com). We would be happy to hear your thoughts! 

[Apache License 2.0](https://github.com/rodrgo/openph/blob/master/LICENSE)

