import os
import seq2simulate

class Platform:

    def __init__(
        self,
        coverage,
        mean_read_length,
        prevalence_error,
        profile=None
    ):
        self.coverage = coverage
        self.profile = profile
        self.mean_read_length = mean_read_length
        self.prevalence_error = prevalence_error

package_dir = seq2simulate.__path__[0]
roche_profile_dir = os.path.join(package_dir, 'profiles/roche')
ion_profile_dir = os.path.join(package_dir, 'profiles/ion')
pacbio_ccs_profile = os.path.join(package_dir, 'profiles/model_qc_ccs')
pacbio_clr_profile = os.path.join(package_dir, 'profiles/model_qc_clr')

illumina = Platform(10000, 250, 0.01)
roche = Platform(3000, 320, 0.03, profile=roche_profile_dir)
ion = Platform(5000, 320, 0.03, profile=ion_profile_dir)
pacbio_ccs = Platform(10000, 250, 0.02, profile=pacbio_ccs_profile)
pacbio_clr = Platform(10000, 250, 0.10, profile=pacbio_clr_profile)

platforms = {
    'roche': roche,
    'illumina': illumina,
    'ion': ion,
    'pacbio-ccs': pacbio_ccs,
    'pacbio-clr': pacbio_clr
}

platform_names = { v: k for k, v in platforms.iteritems() }