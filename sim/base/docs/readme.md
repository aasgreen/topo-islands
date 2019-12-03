# Landau-Ginzburg Simulations of 2D XY-Model with Variable Radius
Contributors:
1. Adam Green (github.com/aasgreen)

The document is organized in the following way. 
1. First, the experimental background is briefly discussed to give some context.
2. The goals of the project is discussed
3. A very basic guide to running the simulations and programs
4. More in-depth guide to each program

## Background

In experimental systems of freely-suspended films, we observed the creation of topological defect pairs that lived far longer than usual. These defect pairs were created by the nucleation of islands in the film during a quench. The islands would grow, creating the defects, and then shrink. Suprinsingly, when the islands disappeared, the defects lived on, proving to be unusually stable.

## Goals

This project was designed to simulate the above system using fortran code previously written in the course of my PhD defense.


## Up and Running
I will be providing example code of what I run on my system (ubuntu on windows 2), this may or may not work for you. 

This will outline the very basic steps to running your first simulation with this code. Because we are using docker, you first have to download/install docker for your respective operating system. (Docker is a 'container', where all the dependancies required to run the programs are bundled, so the users don't have to worry about this). 
### Setup Environment

```
sudo service docker start
```
Then, you will need to pull the latest docker image from my docker file, and run it. This functionality is contained within the "docker_run.sh" script, so run that as follows

```
chmod +x docker_run.sh
./docker_run.sh
```

This should also launch you into a docker container under the directory /work, which is a shared directory between the docker container and the root directory of this project.

### Running your first simulation

The workhorse of this simulation is the program lg.f90, which will setup a grid of spins, and put a hard circularly boundary somewhere in the middle of it. It will then evolve the program using a first-order discretized ODE solver (the euler method), which steps the system through it's equations of motion. You will need to pass the following parameters:

------------------------
Parameters| Description|
------------------------
elastic constant (k) | The interaction energy between neighboring spins[1]|
temperature (beta) | This actually sets the inverse temperature 1/T of the system |
external field (mu) | If you want to apply an external field that acts to globally align the spins, you can do that here |
gridSize | The length of the grid |
endTime | The length of time you want the simulation to evolve for |
seed    | For some operating systems (not unix based), fortran needs an external seed to start it's number generator, so you can provide one here. If you are running this in docker, it won't be used|
----------------------------------------------------------------------

To use the program, first compile it with gfortran:

```
gfortran -03 -o defect.o simulations/lg.f90
```

and then call the program using the parameters you want:

```
./defect.o 1.00 2.00 0.00 100 100 101483
```

So, your current working directory should be overwhelmed with two types of files:

defect$number.dat and out$number.dat

The $number is the time since the beginning of the simulation. The files named defect contain our best guess as to where the defects are located (which is done using an alogorithmn). They are labelled with their respective winding number. The out file is the actual spin grid, containing the spin in radians of each grid-site. 

In practice, I've built a handler script that will put all these data files in a tagged folder, so you don't have to mess up your nice clean initial directory. So if you delete all those data files:

```
rm *.dat
```

And instead run:
```
bash simple_run.sh
```

You should see a new folder, called dataFolder. Inside that should be a folder with the prefix run, and then a long number. This number is the time stamp of your run date, to help keep different simulation runs organized. Inside that folder you should have a few things. The local files should be copied over, so you can rerun this particular simulation run again if you need to. You should have a folder with all the parameters labeled in its name, containing all the raw data.

Unpacking what simple_run.sh does seems a little overwhelming, but its main purpose is to handle directory creation and organization for you. It also calls a very simple python script, which can handle generating lists of parameters (say you wanted to look at different temperatures, you simply make a list of all the temperatures that you want (here I use the np.linspace() function). It will combine that list of temperatues with all the other lists of parameters you want,
creating a text file where each line is parameter set to run the program over. It creates directories, enters them, and then loops over that parameter file, feeding each line into the lg.f90 program.

It has a bunch of lines commented out. A lot of these commented out lines handle more complex functionality, like creating videos and images of the simulation.


## Complex Run

This functionality is based around the bash script call_defect.sh. It has the same basic procudure as the simple_run.sh setup, but it will additionally render movies and images of the simulation as schlieren textures (this is what the spins would look like when viewed under cross polarizers).

The additional moving parts are conceptually simple, even if that simplicity is hidden in the code. Basically, taking a spin and converting it into an image intensity is done through the schlieren map: f(x) = sin^4(x), where x is the angle of the spin. However, dealing with images and annotating them appropriately is a huge pain, at least when I wrote this, which needed a lot of boilerplate.











