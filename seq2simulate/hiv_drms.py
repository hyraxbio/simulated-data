import csv
import os

import drm
import seq2simulate

package_dir = seq2simulate.__path__[0]
drm_tsv = os.path.join(package_dir, 'all_drms.tsv')
drms = []

with open(drm_tsv, 'rb') as csvfile:
	drm_strings = csv.reader(csvfile, delimiter='\t')
	drms = set([drm.Drm(d[1], drm.loci[d[2]]) for d in drm_strings])

