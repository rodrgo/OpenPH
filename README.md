# OpenPH: Persistent Homology computation on GPUs

OpenPH is a library for fast persistent homology computation on GPUs. OpenPH provides a `CUDA-C` implementation of `pms`, which is a provably convergent parallel algorithm for boundary matrix reduction tailored for GPUs, see [Arxiv](https://arxiv.org/abs/1708.04710). Additionally, it includes vanilla implementations of other classical reduction algorithms like `standard`, `twist` and `ph-row`.

## Installation

### OpenPH

The following programs are required to build OpenPH: `make`, `python`, `nvcc`, `wget`, `matlab`. At this moment, `matlab` is the only interface for OpenPH. However, `python` and `julia` interfaces will be available in release `1.0.0`.

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

The following libraries will be installed as the are needed to run the numerical experiments. Their dependencies are listed below each library.

* [Javaplex](https://github.com/appliedtopology/javaplex) - To generate pointclouds, simplicial complexes and boundary matrices.
    * `java`, `ant`
* [Phat](https://bitbucket.org/phat-code/phat) - Benchmarking baseline.
    * `cmake`, `make`
* [DIPHA](https://github.com/DIPHA/dipha) - Benchmarking datasets.
    * `cmake`, `make`, `openmp`

## Contributing

The contributing guidelines will be established in release `1.0.0`. If in the meantime, you want to contribute, please fork the project, create a new branch, and 

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

[MIT License](https://github.com/rodrgo/openph/blob/master/LICENSE)

