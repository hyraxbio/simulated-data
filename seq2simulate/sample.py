import drm
import json

resistance_scores = [(60, 'R'), (15, 'I')]

# There's a problem calling resistance from M184V using Sierra,
# namely that if there's any M184M in the mix, the call goes up to
# resistant.  This only happens in NGS data, so Sierra can't capture
# it and we must treat the case manually here.
fixed_calls = {Drm('M184V', drm.RT) : {'D4T': 'R', 'AZT': 'R'}}

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
		'mutations': [str(d) for d in hiv_drms.drms]
	}

class Sample:
	def __init__(self, name, sequence, prevalence):
		self.name = name
		self.sequence = sequence
		self.prevalence = prevalence

	def dump_csv(self):
		return ','.join(
			self.name,
	        self.sequence.susceptible.id,
	        self.sequence.resistant.id,
	        self.prevalence,
	        str(self.sequence.years_infected),
	        str(self.sequence.pcr_error),
	        str(self.sequence.human_error),
	        str(self.sequence.env_error),
	        str(self.sequence.remove_rt)
        ]) + '\n'

	def acceptable_error(drm):
		# insertion mutations can be up to 10% below their prevalence.
		if (drm.mutation == 'i'):
			return [10, 0]
		else:
			return [0, 0]

	def encode_calls(raw_calls):
		final_calls = {}
		for drug, score in raw_calls:
			found = False
			for drm, drugs_affected in fixed_calls.iteritems():
				if drm in self.sequence.drms \
				and drug in drugs_affected.keys():
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
    	
		mutations = { str(drm) : {
				'prevalence': self.prevalence,
				'acceptable_error': acceptable_error(drm)
			}
		} for drm in self.drms 

    	calls = encode_calls(sierra_wrapper.get_calls(self.sequence.drms))

    	return {'mutations': mutations, 'calls': calls}
