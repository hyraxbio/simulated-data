import csv
import os
import pickle

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

hiv_drms_cache_file = ".hiv_drm_cache"

get_drms_from_service = True
if os.path.isfile(hiv_drms_cache_file):
    print "DRM cache found, checking version."
    with open(hiv_drms_cache_file, 'r') as f:
        cached_version = f.readline().strip()
        ws_version = sierra.get_version()
        print "DRM cache is version", cached_version + ", web service is version", ws_version
        # if there hasn't been a version bump, we can use the cached drms
        if cached_version == ws_version:
            print "Importing from cache."
            get_drms_from_service = False
            drms = [drm.Drm(line) for line in f]

if get_drms_from_service:
    print "Creating new cache from web service."

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

    with open(hiv_drms_cache_file, 'w') as f:
        f.write(sierra.get_version() + '\n')
        for d in drms:
            f.write(str(d) + '\n')