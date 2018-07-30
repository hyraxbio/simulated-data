import os
try:
    from setuptools import setup
    from setuptools.command.install import install
except ImportError:
    from distutils.core import setup
    from distutils.command.install import install



class CheckPresenceOfExecutables(install):
    """
    Customized setuptools install command - checks for presence of
    necessary executables in path.
    """

    # needed_executables = {
    #     'EvolveAGene4': 'https://popmodels.cancercontrol.cancer.gov/' \
    #                     'gsr/packages/evolveagene/',
    #     'java': 'http://www.oracle.com',
    #     'art_454': 'http://www.niehs.nih.gov/research/resources/software/biostatistics/art/',
    #     'art_illumina': 'http://www.niehs.nih.gov/research/resources/software/biostatistics/art/',
    #     'pbsim': 'https://code.google.com/archive/p/pbsim/'
    # }

    needed_executables = {} 

    class ExecutableException(Exception):
        pass

    def which(self, program):
        """
        Check if an executable exists.

        Args:
            program: An executable name or path.

        Returns:
            The full path to the program on success, else None.
        """


        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

        return None


    def run(self):
        for executable in self.needed_executables.keys():
            if self.which(executable) is None:
                raise self.ExecutableException(
                    "Executable " + executable + " is not present. " \
                    "Please download and install it from "
                    + self.needed_executables[executable]
                )
        install.run(self)

setup(
    name='simulate',
    description='Simulated data generator for HIV drug resistance testing',
    author='Imogen Wright',
    author_email='imogen@hyraxbio.co.za',
    version='1.14.0',
    packages=['seq2simulate', 'seq2simulate/revolver', 'seq2simulate/hypermutation'],
    scripts=[
        'bin/all-scores',
        'bin/simulate', 
        'bin/package-simulation', 
        'bin/sierra',
        'bin/merge-samples',
        'bin/make-proviral-data',
    ],
    package_data={'seq2simulate': [
        'all_drms.tsv',
        'contamination/*',
        'profiles/model*',
        'profiles/ion/*',
        'profiles/roche/*',
        'stanford/*.jar',
        'stanford/run.sh',
        'stanford/sierra.properties',
        'stanford/filesNeededBySierra/*.xml',
        'stanford/filesNeededBySierra/*.mat',
        'stanford/filesNeededBySierra/*.txt',
        'stanford/filesNeededBySierra/XMLVirusDefs/*'
    ]},
    include_package_data=True,
    cmdclass={
        'install': CheckPresenceOfExecutables
    },
    install_requires=[
        'appdirs',
        'BioPython',
        'click',
        'pysam',
        'soappy',
        'multiprocessing',
        'wstools',
        'sierrapy',
        'numpy',
        'scipy',
    ],
    setup_requires=[
        'appdirs',
        'BioPython',
        'click',
        'pysam',
        'soappy',
        'multiprocessing',
        'wstools',
        'sierrapy',
        'numpy',
    ]
)

