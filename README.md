# Parallel reduction for persistent homology

This library implements a proof of concept for a parallel algorithm for persistent homology.

## Getting Started

See below for prerequisites and installation instructions.

### Prerequisites

This library requires [MATLAB](http://www.mathworks.com/products/matlab.html) and a particular fork of the [javaplex](https://github.com/appliedtopology/javaplex) library.

This fork can be cloned by

```
git clone https://github.com/rodrgo/javaplex.git
```

### Installing

To install, go to the `javaplex` project folder and compile.

```
cd PATH_TO_JAVAPLEX
ant compile jar run 
```

Then, go to the `tda` project folder and edit the `JAVAPLEX_DIR` variable in the `load_javaplex.m`

```
JAVAPLEX_DIR = 'PATH_TO_JAVAPLEX_DIRECTORY'
```

You are good to go!

## Running the tests

To check that everything is running correctly,

```
cd tests
matlab -nodisplay -nodesktop -r "run tests.m"
```

## Built With

* [Matlab](http://www.mathworks.com/products/matlab.html) - The language used 
* [Javaplex](https://github.com/appliedtopology/javaplex) - To generate of boundary boundary matrices

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Rodrigo Mendoza-Smith** - [personal webpage](people.maths.ox.ac.uk/mendozasmith/)
* **Jared Tanner** - [personal webpage](people.maths.ox.ac.uk/tanner/)

## Citing

TBD

## License

TBD  - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Javaplex library
* Vidit Nanda

