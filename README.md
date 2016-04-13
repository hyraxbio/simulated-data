# simulated-data

## About

Create simulated HIV drug resistance testing data.

Turns a set of susceptible and a set of resistant HIV pol sequences into simulated sample data, as input to a drug resistance testing platform.  Produces a folder containing hashed filenames, and a manifest containing a legend to decode the files.

Example input can be found in the ```data``` folder.

Uses the EvolveAGene4 evolution simulator, the ART sequence simulator and the sierra HIV drug resistance database.


## Install

From root:

### Global
```
python setup.py build
sudo python setup.py install
```

### Local
```
virtualenv env
env/bin/pip install .
```
The install script may suggest a number of external software packages which must be in your build system for seq2simulate to work.  These include EvolveAGene4 and art_454/art_illumina.

* EvolveAGene4: https://popmodels.cancercontrol.cancer.gov/gsr/packages/evolveagene/
* ART: http://www.niehs.nih.gov/research/resources/software/biostatistics/art/

Java is also required, and can be installed via a package manager or from [Oracle](http://www.oracle.com).

## Basic Usage

Create a Roche/454 simulation from Seq6 at 1% prevalence:

```
simulate --prevalence 1 --susceptible data/split/Seq6_Sus --resistant data/split/Seq6_Res --out temp_roche --platform roche
```
Create an Illumina simulation from all sequences, with randomized errors, at all prevalences:

```
env/bin/simulate --susceptible data/SusceptibleSeqs_all.fasta --resistant data/ResistantSeqs_all.fasta --out temp_out --platform illumina --paired-end --randomize
```


# Options:
```
  --susceptible PATH              Path to a FASTA file of susceptible
                                  sequence(s).

  --resistant PATH                Path to a FASTA file of resistant
                                  sequence(s).

  --out PATH                      Path to the desired output folder (will be
                                  created).

  --platform [illumina|ion|roche]
                                  Sequencing platform to simulate.

  --working-dir PATH              Returns a specific working directory
                                  location. Defaults to a randomly generated
                                  uuid in an appropriate data folder.

  --prevalence INTEGER            Select a prevalence to run at.  Default is
                                  to generate 10prevalences: 0.5, 1, 2, 5, 10,
                                  15, 20, 30, 50, 100

  --diversity-only                Leaves the working directory at a point
                                  where diversity, but not sequencing error,
                                  has been simulated.

  --error-only                    When used with --working-dir, generates
                                  sequencing error from the given working
                                  directory, which should already have been
                                  created using --evolve-only.

  --prevalence-only               When used with --working-dir, creates the
                                  final output from already --error-only
                                  simulated reads in the given working
                                  directory

  --pcr-error                     Include a false PCR error in this
                                  simulation.

  --env-error                     Include some DNA from HIV ENV in this
                                  simulation.

  --human-error                   Include some human DNA in this simulation.
  --remove-rt                     Remove the piece of RT covering K103N
                                  entirely.

  --randomize                     For multiple sequences.  Shuffle and
                                  randomly assign error types.

  --paired-end                    Paired-end sequencing.  Only valid with
                                  illumina data.

  --help                          Show this message and exit.
```
