
import appdirs
import json
import ntpath
from multiprocessing.pool import Pool
import multiprocessing
import os
import parser
import pickle
import Queue
import random
import time
import traceback

import art
import diversity
import prevalence
import sample
import sequencing_error

evolved_filename = 'evolved_files.pkl'
error_filename = 'error_files.pkl'
final_csv_filename = 'manifest.csv'
final_json_filename = 'manifest.json'
max_years_infected = 10

csv_header_rows = [
    'Hashed_filename',
    'Susceptible',
    'Resistant',
    'DRM_prevalence',
    'Years_infected',
    'PCR_error',
    'Human_DNA',
    'ENV_DNA',
    'RT_Removed'
]

platforms = {
    'roche': art.roche,
    'illumina': art.illumina,
    'ion': art.ion,
    'pacbio-ccs': art.pacbio_ccs,
    'pacbio-clr': art.pacbio_clr
}

DEFAULT_PREVALENCES = [0.005, 0.01, 0.02, 0.05, 0.10, 
                        0.15, 0.20, 0.30, 0.50, 1]

MINIMUM_FOR_SHUFFLE = 4

# Multiprocessing in Python 2 has buggy error reporting.  This
# gives us a stack trace
# Thanks to http://stackoverflow.com/questions/6728236/exception-thrown-in-multiprocessing-pool-not-detected
class LogExceptions(object):
    """
    Class that logs exceptions from multiprocessing threads.
    """
    def __init__(self, callable):
        self.__callable = callable
        return

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)

        except Exception as e:
            # Here we add some debugging help. If multiprocessing's
            # debugging is on, it will arrange to log the traceback
            multiprocessing.get_logger().error(traceback.format_exc())
            # Re-raise the original exception so the Pool worker can
            # clean up
            raise

        # It was fine, give a normal answer
        return result
    pass

class LoggingPool(Pool):
    """
    Extend pool and overload apply_async
    """
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), 
            args, kwds, callback)

def run_diversity_thread(evolved_queue, sequence, working_dir):

    """
    Calculates diversity for a single sequence

    Args:
        evolved_data: The list of calculated diversity data.
        sequence: The sequence to simulate.
        working_dir: Put temp files here.
    """
    
    files = diversity.simulate(sequence, working_dir)
    files["sequence"] = sequence
    evolved_queue.put(files)

def run_diversity(
    raw_susceptible, raw_resistant,
    pcr_error, env_error, human_error, remove_rt,
    randomize,
    working_dir, 

):
    """
    Simulate diversity from sets of susceptible and resistant sequences.

    Args:
        raw_susceptible: A filename containing susceptible sequences.
        raw_resistant: A filename containing their equivalent resistant
            sequences.
        pcr_error, env_error, human_error: Force these error types.
        randomize: Assign random error types.
        working_dir: Put temp files here.

    Returns:
        True on completion.
    """

    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    sequences = parser.parse_sequences(raw_susceptible, raw_resistant)

    if randomize:
        if len(sequences) < MINIMUM_FOR_SHUFFLE:
            raise ValueError(
                "Not enough sequences for random assignment. "\
                "Need " + str(MINIMUM_FOR_SHUFFLE) + "."
            )
        if pcr_error or human_error or env_error:
            raise ValueError(
                "Cannot assign both random and specific errors."
            )
        print "Shuffling and randomly assigning error types."

        sequencing_error.shuffle_and_assign_error_types(sequences)

    
    threads = []
    evolved_queue = multiprocessing.Manager().Queue()

    multiprocessing.log_to_stderr()
    p = LoggingPool(len(sequences))

    for sequence in sequences:

        # each patient has been infected for a number of years with the
        # original infection being the given sequence
        sequence.years_infected = random.randint(1, max_years_infected)

        if pcr_error:
            sequence.pcr_error = True
        if env_error:
            sequence.env_error = True
        if human_error:
            sequence.env_error = True
        if remove_rt:
            sequence.remove_rt = True

        p.apply_async(run_diversity_thread, [
            evolved_queue, sequence, working_dir
        ])
    
    p.close()
    p.join()

    evolved_data = []
    while not evolved_queue.empty():
        evolved_data.append(evolved_queue.get())

    with open(os.path.join(working_dir, evolved_filename), 'wb') as handle:
        pickle.dump(evolved_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return True


def run_error_thread(error_queue, result, platform, working_dir, 
    pcr_error, env_error, human_error, paired_end):

    """
    Simulate for a single sequence

    Args:
        error_data: Where to put results.
        result: The current diversity result to process.
        platform: One of "roche", "illumina" or "ion".
        working_dir: The folder in which to place temporary files.
        pcr_error: Should we include a PCR error?
        env_error: Should we include an ENV error?
        human_error: Should we include human DNA?

    Returns:
        True on completion.
    """

    error_result = {'platform' : platform}

    error_result['sequence'] = result['sequence']
    if pcr_error:
        result['sequence'].pcr_error = True
    if env_error:
        result['sequence'].env_error = True
    if human_error:
        result['sequence'].human_error = True

    error_result['susceptible_fq'], error_result['susceptible_sam'] =\
        sequencing_error.simulate(
            result['sequence'], result['susceptible'], 
            platforms[platform], paired_end, working_dir
    )
    error_result['resistant_fq'], error_result['resistant_sam'] =\
        sequencing_error.simulate(
            result['sequence'], result['resistant'], 
            platforms[platform], paired_end, working_dir
    )

    error_queue.put(error_result)

def run_error(platform, working_dir, pcr_error, env_error, 
    human_error, paired_end):
    """
    Simulate sequencing error from the output generated by run_diversity.

    Args:
        platform: string One of "roche", "illumina" or "ion".
        working_dir: path The folder in which to place temporary files.
        pcr_error: bool Should we include a PCR error
        env_error: bool Should we include an ENV error
        human_error: bool Should we include human DNA
        paired_end: bool Are we simulating paired_end data?

    Returns:
        True on completion.
    """

    print "Simulating reads from sequence sets."

    error_data = {}
    # properties of the simulation
    error_data['paired_end'] = paired_end
    error_data['platform'] = platform

    with open(os.path.join(working_dir, evolved_filename), 'rb') as handle:

        evolved_data = pickle.load(handle)

        multiprocessing.log_to_stderr()
        p = LoggingPool(len(evolved_data))
        file_queue = multiprocessing.Manager().Queue()
       
        for result in evolved_data:

            p.apply_async(run_error_thread, [
                file_queue, result, platform, working_dir, 
                pcr_error, env_error, human_error, paired_end
            ])

        p.close()
        p.join()

        error_data['files'] = []
        while not file_queue.empty():
            error_data['files'].append(file_queue.get())

    with open(os.path.join(working_dir, error_filename), 'wb') as handle:
        pickle.dump(error_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return True


def run_prevalence_thread(manifest_queue, platform, paired_end,
    result, working_dir, out_dir, produce_prevalence):
    """
    Single-threaded generator of pools of reads.

    Args:
        manifest_queue: The parallel queue onto which to push manifests.
        platform: One of "roche", "illumina" or "ion".
        paired_end: Are we simulating paired_end data?
        result: The sequencing error result from which to simulate prevs.
        working_dir: The folder in which to place temporary files.
        out_dir: The final output directory to place the file.
        prevalence: The prevalence at which to run (default is all)
        
    Returns:
        True on completion.
    """

    print "Simulating reads from sequence sets."

    prevalences = DEFAULT_PREVALENCES
    if produce_prevalence is not None:
        prevalences = [produce_prevalence]

    for required_prevalence in prevalences:

        hashed_filename = prevalence.produce_prevalence(
            required_prevalence,
            platforms[platform],
            paired_end,
            result['susceptible_fq'], result['susceptible_sam'],
            result['resistant_fq'], result['resistant_sam'],
            result['sequence'],
            working_dir,
            out_dir
        )

        manifest_queue.put(sample.Sample(
            ntpath.basename(hashed_filename),
            result['sequence'],
            required_prevalence
        ))
    

def run_prevalence(out_dir, remove_rt, working_dir, produce_prevalence):
    """
    Create exact prevalences over a list of DRMs provided by a sequence
    object when provided with simulated resistant and susceptible files.

    Args:
        out_dir: The final folder into which to place the anonymized files.
        remove_rt: Are we removing a piece of the RT gene?
        working_dir: The folder in which to place temporary files.
        produce_prevalence: what prevalence to produce

    Returns:
        True on completion.
    """

    print "Building exact prevalences."

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    csv_rows = []
    threads = []

    with open(os.path.join(working_dir, error_filename), 'rb') as handle:
        error_data = pickle.load(handle)
        platform = error_data['platform']
        paired_end = error_data['paired_end']

        manifest_queue = multiprocessing.Manager().Queue()
        multiprocessing.log_to_stderr()
        p = LoggingPool(len(error_data['files']))
    
        for result in error_data['files']:

            if remove_rt:
                result['sequence'].remove_rt = True

            process_result = p.apply_async(run_prevalence_thread, [
                manifest_queue, platform, paired_end,
                    result, working_dir, out_dir, produce_prevalence
            ])

        p.close()
        p.join()

        with open(os.path.join(out_dir, 
            final_csv_filename), 'w') as csv_handle,\
            open(os.path.join(out_dir, 
            final_json_filename), 'w') as json_handle:
            csv_handle.write(','.join(csv_header_rows) + '\n')
            test = sample.header(platforms[error_data['platform']])
            test['samples'] = []
            while not manifest_queue.empty():
                s = manifest_queue.get()
                test['samples'].append(s.encode())
                csv_handle.write(s.dump_csv())
            json_handle.write(json.dumps(test, indent=2))
                
    return True