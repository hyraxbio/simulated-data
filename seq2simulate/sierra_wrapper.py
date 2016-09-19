# WSDL for the service: http://hivdb.stanford.edu/DR/schema/sierra.xsd

import os
import uuid
import pkgutil
import Bio.SeqIO
import shutil
import SOAPpy
import subprocess
import xml.etree.ElementTree as ET

import drm
import seq2simulate


key = 'B0VQ-C341-5FPU-LO2P' 
endpoint = 'http://db-webservices.stanford.edu:5440/' \
        'axis/services/StanfordAlgorithm?wsdl'

executable = 'run.sh'

server = SOAPpy.SOAPProxy(endpoint)
report_version = 2 # a clearer report version with a defined schema 

# slightly different from definition in drm.py
mutation_dict = {
    'PR' : drm.PR, 
    'RT' : drm.RT,
    'IN' : drm.INI
}

combination_mutations = ['E44A', 'M230I', 'V118I']

package_dir = seq2simulate.__path__[0]
sierra_dir = os.path.join(package_dir, 'stanford')

class SierraError(Exception):
    pass

def call_sierra(sequence):
    """
    Call local sierra JAVA instance.

    Args:
        sequence: A Biopython sequence object.

    Returns:
        Flat text produced by Sierra.
    """

    unique = str(uuid.uuid4())

    temp_seq_filename = os.path.join("/tmp", 
        unique + '_sierra_temp_file.fasta')
    temp_out_filename = os.path.join("/tmp",
        unique + '_sierra_report_file.xml')

    with open(temp_seq_filename, 'w') as handle:
        Bio.SeqIO.write([sequence], handle, 'fasta')

    arg_list = ['/bin/bash', os.path.join(sierra_dir, executable), 
        temp_seq_filename,
        temp_out_filename]

    subprocess.call(arg_list)

    with open(temp_out_filename, 'r') as handle:
        sierra_result = ''.join(handle.readlines())

    os.unlink(temp_seq_filename)
    os.unlink(temp_out_filename)

    return sierra_result

def call_sierra_with_retries(sequence):
    """
    Call sierra until success or N retries and return document root of XML.

    Args:
        sequence: A Biopython sequence object.

    Returns:
        Sierra XML documentRoot ET object.
    """
    retries = 5
    while retries > 0:
        result_string = call_sierra(sequence)
        root = ET.fromstring(result_string)
        if root.find('result') is None \
        or root.find('result').find('success') is None \
        or root.find('result').find('success').text == 'false':
            retries -= 1
            print 'Sequence could not be processed by server.  Retrying.'
            continue
        break
    if retries == 0:
        raise SierraError('Sequence could not be processed by server')

    return root

def get_calls(sequence):
    """
    Use HIVdb to find drug calls from DRMs

    Args:
        sequence: a Biopython sequence object.

    Returns:
        A dictionary of {"drug_name": score}
    """

    calls = {}
    root = call_sierra_with_retries(sequence)

    for drug_score in root.find('result').findall('drugScore'):
        drug_code = drug_score.find('drugCode').text
        score = drug_score.find('score').text
        if drug_code is not None and score is not None:
            calls[drug_code] = int(score)

    return calls

def get_drms(sequence):
    """
    Use HIVdb to find which DRMs are present in a sequence 

    Args:
        sequence: a BioPython sequence object.

    Returns:
        A list of Drm objects
    """

    drms = []
    
    root = call_sierra_with_retries(sequence)
    for gene_data in root.iter('geneData'):
        gene = gene_data.find('gene')
        for mutation in gene_data.iter('mutation'):
            add = False
            mut_string = mutation.find('mutationString')
            # important mutations that only score in combination
            if gene.text == 'RT' and mut_string.text in combination_mutations:
                add = True
            elif mut_string is not None and gene is not None \
              and gene.text in mutation_dict:
                for drug_score in root.iter('drugScore'):
                    for partial_score in drug_score.iter('partialScore'):
                        mut = partial_score.find('mutation')
                        #print mut.text
                        if mut is not None and mut.text == mut_string.text:
                            score = partial_score.find('score')
                            if score is not None:
                                if int(score.text) != 0:
                                    add = True

            if add:

                mut_text = mut_string.text
                mut_text = mut_text.replace(' insertion', 'i')
                mut_text = mut_text.replace(' deletion', 'd')
                drms.append(drm.Drm(mut_text, mutation_dict[gene.text]))
    return drms


