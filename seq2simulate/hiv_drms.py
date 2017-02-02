import csv
import os

import drm
import seq2simulate
import sierra_ws as sierra

amino_acids = "ARNDCQEGHILKMFPSTWYV"
# amino acid pr start position + 1
start = 57
# end of pol
end = 1003

# doesn't matter what this is, it's just a placeholder for "wild type"
placeholder = 'A'

potential_drms = []
results = set()
for a in amino_acids:
    drms_with_amino = [
        drm.Drm(placeholder + str(i) + a) for i in range(start, end)]
    potential_drms.extend(drms_with_amino)
    # first, we try each amino individually to avoid max rules
    results = results.union(set(sierra.is_drm(drms_with_amino)))

# then we check them all together to find drms that score only in combination
results = results.union(set(sierra.is_drm(potential_drms)))

# then we do indel mutations separately
for token in 'id':
    results = results.union(set(sierra.is_drm([
        drm.Drm(placeholder + str(i) + token) for i in range(start, end)])))

drms = results

