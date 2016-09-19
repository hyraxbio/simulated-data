import json

import art
import drm
import hiv_drms
import prevalence
import sierra_wrapper

resistance_scores = [(60, 'R'), (15, 'I')]

# There's a problem calling resistance from M184V using Sierra,
# namely that if there's any M184M in the mix, the call goes up to
# resistant.  This only happens in NGS data, so Sierra can't capture
# it and we must treat the case manually here.
fixed_calls = { drm.Drm('M184V', drm.RT) : {'D4T': [50,'R'], 'AZT': [50,'R']}}

# If we remove everything in RT barring K65, we should get low coverage
# calls for all NRTIs and NNRTIs unless we see K65R
low_coverage_drugs = ['ABC', 'DDI', 'FTC', '3TC', 'D4T', 
                      'TDF', 'AZT', 'EFV', 'ETR', 'NVP', 'RPV']
call_if_resistant = ['TDF', 'DDI']


def header(platform):
    result = {}
    result['platform'] = {
        'name': art.platform_names[platform],
        'coverage_depth': platform.coverage,
        'error_bars': [platform.prevalence_error, 
                             platform.prevalence_error]
    }
    result['pathogen'] = {
        'name': 'HIV',
        'mutations': ["[pol] " + str(d) for d in hiv_drms.drms]
    }
    return result

class Sample:
    def __init__(self, name, sequence, prevalence):
        self.name = name
        self.sequence = sequence
        self.prevalence = prevalence

    def dump_csv(self):
        return ','.join([
            self.name,
            self.sequence.susceptible.id,
            self.sequence.resistant.id,
            str(self.prevalence),
            str(self.sequence.years_infected),
            str(self.sequence.pcr_error),
            str(self.sequence.human_error),
            str(self.sequence.env_error),
            str(self.sequence.remove_rt)
        ]) + '\n'

    def is_no_coverage(self, drm):
        return (
            self.sequence.remove_rt \
            and drm.nucleotide_pos >= prevalence.rt_no_coverage_start \
            and drm.nucleotide_pos <= prevalence.rt_no_coverage_end
        )

    def required_prevalence(self, drm):
        if self.is_no_coverage(drm):
            return 0
        else:
            return self.prevalence

    def acceptable_error(self, drm):
        # insertion mutations can be up to 10% below their prevalence.
        if (drm.mutation == 'i'):
            return [0.1, 0]
        else:
            return [0, 0]

    def encode_calls(self, raw_calls):
        final_calls = {}
        for drug, score in raw_calls.iteritems():
            found = False
            for drm, drugs_affected in fixed_calls.iteritems():
                if drm in self.sequence.drms and drug in drugs_affected.keys() \
                and score >= drugs_affected[drug][0]:
                    final_calls[drug] = drugs_affected[drug][1]
                    found = True
            if not found:
                for cutoff, call in resistance_scores:
                    if score >= cutoff:
                        final_calls[drug] = call
                        found = True
                        break
            if not found:
                final_calls[drug] = 'S'

            # check for low coverage
            if (
                self.sequence.remove_rt \
                and drug in low_coverage_drugs \
                and not (
                	final_calls[drug] != 'S' and drug in call_if_resistant
                )
            ):
                final_calls[drug] = 'LC'

        return final_calls


    def encode(self):
        
        mutations = [ 
                {'name': "[pol] " + str(drm),
                'prevalence': self.required_prevalence(drm),
                'error_bars': self.acceptable_error(drm)} \
                for drm in self.sequence.drms ]

        calls = self.encode_calls(sierra_wrapper.get_calls(
            self.sequence.resistant))

        notes = {"simulated_from" : self.sequence.resistant.id}
        notes["errors"] = []
        if self.sequence.remove_rt:
        	notes["errors"].append("RT removed")
        if self.sequence.pcr_error:
        	notes["errors"].append("PCR error")
        if self.sequence.human_error:
        	notes["errors"].append("Human DNA added")
        if self.sequence.env_error:
        	notes["errors"].append("ENV DNA added")

        return {
                'name': self.name,
                'mutations': mutations, 
                'calls': calls,
                'notes': notes
                }
