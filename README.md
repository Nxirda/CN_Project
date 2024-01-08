# Poisson 1D Equation Solver 

## Build with docker :

### Docker :

You will need docker installed on your machine : 

- Docker documentation is available [here](https://docs.docker.com/engine/install/)

### Build the project :

```sh
#Clone this repository
$ git clone https://github.com/Nxirda/CN_Project.git

#First build the image, once in the project directory :
$ docker build buildenv -t cn_config

#Then start the container
$ docker run --rm -it -v "$(pwd)":/root/env cn_config
```

Then you can skip to the Usage section below


## Build : Standard way

### Dependencies :

You will need :
- gcc
- gnuplot
- build essentials (for the make command)
- lapack
- blas

```sh
#This project has been tested to work for Ubuntu 20.04
$ sudo apt install gcc 
$ sudo apt install gnuplot
$ sudo apt install build-essential

$ sudo apt install libblas-dev
$ sudo apt install liblapack-dev  
$ sudo apt install liblapacke-dev  
```

### Build :

```sh
$ git clone https://github.com/Nxirda/CN_Project.git
$ cd TP_Poisson_C_students_2022
$ make -j
```

## Usage :

### Project :
- Refer to the README in the TP_Poisson_C_students_2022 directory

### Plots :

The plot.p file takes a RESVEC.dat file in bin directory

```sh
$ gnuplot plot.p
```

## Authors :
[Adrien Henrot](https://github.com/Nxirda)