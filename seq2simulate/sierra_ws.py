import Bio
import json
import os
import sys
import uuid
from gql import gql
from sierrapy import SierraClient, fragments, fastareader

from seq2simulate import drm

endpoint = 'https://hivdb.stanford.edu/graphql'

# slightly different from definition in drm.py
mutation_dict = {
    'PR' : drm.PR, 
    'RT' : drm.RT,
    'IN' : drm.IN
}

def get_version():
    """
    Get the current version of HIVdb
    """
    client = SierraClient(endpoint)
    data = client.execute(
        gql('query { viewer { currentVersion { text } } }')
    )
    return data['viewer']['currentVersion']['text']


def parse_drms_from_scores(response):
    """
    Parse a sierra generated json list of mutations.  
    Works without a known list of DRMs, but cannot return DRMs
    that score only in combination.

    Args:
        response: the deserialized json response.
    
    Returns:
        drms: all drms in the response
    """
    drms = []
    for gene in response["drugResistance"]:
        gene_name = gene["gene"]["name"]
        for score in gene["drugScores"]:
            for partial_score in score["partialScores"]:
                for mutation in partial_score["mutations"]:
                    mut_text = mutation["text"]
                    mut_text = mut_text.replace("Insertion", "i")
                    mut_text = mut_text.replace("Deletion", "d")
                    mut = drm.Drm(mut_text, mutation_dict[gene_name])
                    # if the mutation scores, or we know it's one
                    # of the drms that scores in combination
                    if partial_score["score"] != 0:
                        if mut not in drms:
                            drms.append(mut)
    return drms  

def parse_drms(response, known_drms=[]):
    """
    Parse a sierra generated json list of mutations.
    Filters out only DRMs from a known list of DRMs.
    Cannot ab initio discover DRMs.

    Args:
        response: the deserialized json response.
        known_drms: the list of already-known drms
    
    Returns:
        drms: all drms in the response
    """
    drms = []
    for gene in response["alignedGeneSequences"]:
        gene_name = gene["gene"]["name"]
        # grab mutations 
        for mutation in gene["mutations"]:
            mutated_to = mutation["AAs"]
            # handle insertions
            if len(mutated_to) > 1:
                mutated_to = "i"
            if mutated_to == "-":
                mutated_to = "d"
            mut_text = mutation["consensus"] \
                      + str(mutation["position"]) \
                      + mutated_to

            mut = drm.Drm(mut_text, mutation_dict[gene_name])
            if mut in known_drms and mut not in drms:
                drms.append(mut)
    # print(drms)
    return drms                        

def call_sierra_with_drms(drms):
    """
    Call the sierra web service with DRMs.
    """

    client = SierraClient(endpoint)
    response = client.mutations_analysis(drms, 
        fragments.MUTATIONS_ANALYSIS_DEFAULT)
    return response

def is_drm(potential_drms):
    """
    Check which of a list of amino acids is a DRM.

    Args: 
        potential_drms: a list of Drm objects.
    Returns:
        drms: a list of which ones were actually DRMs.
    """
    request_drms = [
        drm.locus_names[d.locus] + ":" \
        + d.relative_str()[1:].replace("i", "ins").replace("d", "del")
         for d in potential_drms
    ]
    response = call_sierra_with_drms(request_drms)
    #json.dump(response, sys.stdout, indent=2)
    drms = []
    if len(response) == 0 or "drugResistance" not in response:
        raise ValueError("Not HIV DRMs.")

    return parse_drms_from_scores(response)

def call_sierra_with_sequence(sequence):
    """
    Call the sierra web service.

    Args:
        sequence: A Biopython sequence object.

    Returns:
        encoded JSON object of response
    """

    unique = str(uuid.uuid4())
    sequence_file = os.path.join("/tmp", 
        unique + '_sierra_temp_file.fasta')

    with open(sequence_file, 'w') as handle:
        Bio.SeqIO.write([sequence], handle, 'fasta')
    client = SierraClient(endpoint)
    query = fragments.SEQUENCE_ANALYSIS_DEFAULT
    with open(sequence_file, 'r') as fp:
        sequences = fastareader.load(fp)
        response = client.sequence_analysis(sequences, query)
        return response

def parse_calls(response):
    """
    Use HIVdb to find drug calls from DRMs
    """
    calls = {}
    for gene in response["drugResistance"]:
        for score in gene["drugScores"]:
            drug_name = score["drug"]["displayAbbr"]
            drug_score = score["score"]
            calls[drug_name.encode("ascii")] = drug_score

    return calls


def get_calls_from_sequence(sequence):
    """
    Use HIVdb to find drug calls from a sequence

    Args:
        sequence: a Biopython sequence object.

    Returns:
        A dictionary of {"drug_name": score}
    """

    response = call_sierra_with_sequence(sequence)
    #json.dump(response, sys.stdout, indent=2)
    if len(response) == 0 or "drugResistance" not in response[0]:
        raise ValueError("Not HIV DNA.")
    
    return parse_calls(response[0])

def get_calls_from_drms(drms):
    """
    Use HIVdb to find drug calls from DRMs

    Args:
        sequence: a Biopython sequence object.

    Returns:
        A dictionary of {"drug_name": score}
    """

    request_drms = [
        drm.locus_names[d.locus] + ":" \
        + d.relative_str()[1:].replace("i", "ins").replace("d", "del")
         for d in drms
    ]
    response = call_sierra_with_drms(request_drms)
    
    if len(response) == 0 or "drugResistance" not in response:
        raise ValueError("Not HIV DRMs.")
    
    return parse_calls(response)

def get_drms(sequence, known_drms=[]):
   
    """
    Use HIVdb to find which DRMs are present in a sequence 

    Args:
        sequence: a BioPython sequence object.

    Returns:
        A list of Drm objects
    """
    response = call_sierra_with_sequence(sequence)
    #json.dump(response, sys.stdout, indent=2)

    if len(response) == 0 or "drugResistance" not in response[0]:
        raise ValueError("Not HIV DNA.")
    return parse_drms(response[0], known_drms=known_drms)


