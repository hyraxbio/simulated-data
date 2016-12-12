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

def call_sierra(sequence):
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


def get_calls(sequence):
    """
    Use HIVdb to find drug calls from DRMs

    Args:
        sequence: a Biopython sequence object.

    Returns:
        A dictionary of {"drug_name": score}
    """

    response = call_sierra(sequence)
    #json.dump(response, sys.stdout, indent=2)
    calls = {}
    if len(response) == 0 or "drugResistance" not in response[0]:
        raise ValueError("Not HIV DNA.")
    for gene in response[0]["drugResistance"]:
        for score in gene["drugScores"]:
            drug_name = score["drug"]["displayAbbr"]
            drug_score = score["score"]
            calls[drug_name.encode("ascii")] = drug_score

    return calls


def get_drms(sequence):
   
    """
    Use HIVdb to find which DRMs are present in a sequence 

    Args:
        sequence: a BioPython sequence object.

    Returns:
        A list of Drm objects
    """
    response = call_sierra(sequence)
    if len(response) == 0 or "drugResistance" not in response[0]:
        raise ValueError("Not HIV DNA.")
    drms = []
    for gene in response[0]["drugResistance"]:
        gene_name = gene["gene"]["name"]
        for score in gene["drugScores"]:
            for partial_score in score["partialScores"]:
                if partial_score["score"] != 0:
                    for mutation in partial_score["mutations"]:
                        mut_text = mutation["text"]
                        mut_text = mut_text.replace("Insertion", "i")
                        mut_text = mut_text.replace("Deletion", "d")
                        mut = drm.Drm(mut_text, mutation_dict[gene_name])
                        if mut not in drms:
                            drms.append(mut)

    return drms