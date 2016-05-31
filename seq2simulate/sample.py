import json
import art
import drm
import hiv_drms
import sierra_wrapper

resistance_scores = [(60, 'R'), (15, 'I')]

# There's a problem calling resistance from M184V using Sierra,
# namely that if there's any M184M in the mix, the call goes up to
# resistant.  This only happens in NGS data, so Sierra can't capture
# it and we must treat the case manually here.
fixed_calls = { drm.Drm('M184V', drm.RT) : {'D4T': 'R', 'AZT': 'R'}}

def header(platform):
    result = {}
    result['platform'] = {
        'name': art.platform_names[platform],
        'coverage_depth': platform.coverage,
        'acceptable_error': [platform.prevalence_error, 
                             platform.prevalence_error]
    }
    result['organism'] = {
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

    def acceptable_error(self, drm):
        # insertion mutations can be up to 10% below their prevalence.
        if (drm.mutation == 'i'):
            return [10, 0]
        else:
            return [0, 0]

    def encode_calls(self, raw_calls):
        final_calls = {}
        for drug, score in raw_calls.iteritems():
            found = False
            for drm, drugs_affected in fixed_calls.iteritems():
                if drm in self.sequence.drms and drug in drugs_affected.keys():
                    final_calls[drug] = drugs_affected[drug]
                    found = True
            if not found:
                for cutoff, call in resistance_scores:
                    if score >= cutoff:
                        final_calls[drug] = call
                        found = True
                        break
            if not found:
                final_calls[drug] = 'S'

        return final_calls


    def encode(self):
        
        mutations = { "[pol] " + str(drm) : {
                'prevalence': self.prevalence,
                'acceptable_error': self.acceptable_error(drm)
            } for drm in self.sequence.drms 
        } 

        calls = self.encode_calls(sierra_wrapper.get_calls(
        	self.sequence.resistant))

        return {'mutations': mutations, 'calls': calls}
