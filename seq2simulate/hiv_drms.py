import csv

import drm

drms = []
with open('all_drms.tsv', 'rb') as csvfile:
	drm_strings = csv.reader(csvfile, delimiter='\t')
	drms = [Drm(d[0], drm.gene_regions[d[2]] for d in drm_strings]

