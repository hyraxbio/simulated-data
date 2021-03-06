import json

import art
import drm
import hiv_drms
import platform as plat
import prevalence
import sierra_ws as sierra

resistance_scores = [(60, 'R'), (15, 'I')]

# There's a problem calling resistance from M184V using Sierra,
# namely that if there's any M184M in the mix, the call goes up to
# resistant.  This only happens in NGS data, so Sierra can't capture
# it and we must treat the case manually here.
fixed_calls_below_100_prevalence = { drm.Drm('M184V', drm.RT) : {'D4T': [50,'R'], 'AZT': [50,'R']}}

# If we remove everything in RT barring K65, we should get low coverage
# calls for all NRTIs and NNRTIs unless we see K65R
low_coverage_drugs = ['ABC', 'DDI', 'FTC', '3TC', 'D4T', 'DOR', 
                      'TDF', 'AZT', 'EFV', 'ETR', 'NVP', 'RPV']
call_if_resistant = ['TDF', 'DDI', 'D4T']


def header(platform, paired_end=False):
    """
    Return the manifest header in sigma format.
    Args:
        platform: the platform Enum from platform
    """
    result = {}

    error = platform.prevalence_error

    # side effect of simulated data's read selection strategy
    # TODO: make this irrelevant
    if paired_end:
        error *= 2

    result['platform'] = {
        'name': plat.platform_names[platform],
        'coverage_depth': platform.coverage,
        'paired_end': paired_end,
        'error_bars': [error, error]
    }
    result['pathogen'] = {
        'name': 'HIV',
        'mutations': [d.locus_str() for d in hiv_drms.drms]
    }
    return result

class Sample:
    """
    A view on the data generated by simulated data.
    """
    def __init__(self, name, sequence, prevalence):
        """
        Constructor.

        Args:
            name: the name of the sample
            sequence: the sequence object associated with the sample
            prevalence: the prevalence at which resistant mutations appear
        """
        self.name = name
        self.sequence = sequence
        self.prevalence = prevalence

    def dump_csv(self):
        """
        Produce a csv manifest.
        """
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
        """
        Returns logic to figure out if the DRM should be covered.
        Right now, DRMs are always covered unless we manually remove
        a piece of RT.
        Args:
            drm: The drug resistance mutation in question
        Returns:
            True if the DRM should not be covered.
        """
        return (
            self.sequence.remove_rt \
            and drm.nucleotide_pos >= prevalence.rt_no_coverage_start \
            and drm.nucleotide_pos <= prevalence.rt_no_coverage_end
        )

    def required_prevalence(self, drm):
        """
        Returns the required prevalence for a DRM.
        Args:
            drm: The drug resistance mutation in question.
        Returns:
            The prevalence of the DRM (usually the prevalence of the sample).
        """
        if self.is_no_coverage(drm):
            return 0
        else:
            return self.prevalence


    def acceptable_error(self, drm, platform):
        """
        Returns the error bars for a particular mutation.
        Args:
            drm: The drug resistance mutation in question.
            platform: the sequencing platform of interest
        Returns:
            A 2 element list, with the lower and upper error bounds
            for the mutation.  Currently this is always zero (and thus superseded
            by the platform error bars) unless the mutation is an indel.
        """

        # insertion mutations can be up to relative 10% below 
        # their prevalence.
        if drm.mutation == 'i':
            return [0.1 * float(self.required_prevalence(drm)), 0]
        # this particular linkage is impossible to call as accurately
        # as many others.
        elif str(drm) == "D222N" and ("T224i" in \
            [str(s) for s in self.sequence.drms]):
            return [0.05 * float(self.required_prevalence(drm)), 0]

        # we're more lenient when the susceptible/resistant ratio
        # approaches 50/50
        elif self.required_prevalence(drm) >= 0.20 and \
             self.required_prevalence(drm) <= 0.80:
            return [
                platform.prevalence_error * 2,
                platform.prevalence_error * 2
            ]

        else:
            return [0, 0]

    def encode_calls(self, raw_calls):
        """
        Encode a sierra resistance score into a call.
        Args:
            raw_calls: The calls from sierra.
        Returns:
            The final calls encoded is a dictionary of:
                drug: {'S' | 'I' | 'R' | 'LC'}
        """
        final_calls = {}
        for drug, score in raw_calls.iteritems():
            found = False

            # some calls need to be artificially doctored at mixed prevalence
            for drm, drugs_affected in fixed_calls_below_100_prevalence.iteritems():
                if drm in self.sequence.drms and \
                   self.required_prevalence(drm) < 1.0 and \
                   drug in drugs_affected.keys() and \
                   score >= drugs_affected[drug][0]:
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
                    final_calls[drug] == 'R' and drug in call_if_resistant
                )
            ):
                final_calls[drug] = 'LC'

        return final_calls


    def make_susceptible(self, call):
        """
        Convert a call to susceptible, unless it's low coverage
        Args:
            call: the call to convert
        Returns:
            The converted call.
        """
        if call != 'LC':
            return 'S'
        return 'LC'


    def encode(self, platform):
        """
        Encode a sample as a dictionary that can be easily JSON-serialized.
        Args:
            platform: the platform for which to set the prevalence range.

        Returns:
            A dictionary.
        """

        if prevalence.range_unusable(platform, self.prevalence):
            raise ValueError("Can't encode this sample for given platform - prevalence ambiguous.")

        mutations = [ 
                {'name': drm.locus_str(),
                'prevalence': self.required_prevalence(drm),
                'error_bars': self.acceptable_error(drm, platform)} \
                for drm in self.sequence.drms ]

        # If the DRM prevalence is low enough that we can test for 
        # "susceptibility"
        if prevalence.range_susceptible(platform, self.prevalence):
            calls = self.encode_calls(sierra.get_calls_from_sequence(
                self.sequence.susceptible))

            prevalence_range = [
                platform.prevalence_error, 
                1
            ]
        # If the DRM prevalence is high enough that we can test for
        # true resistance
        else:
            calls = self.encode_calls(sierra.get_calls_from_sequence(
            self.sequence.resistant))
            
            # The resistant prevalence range runs from the platform error rate to
            # true prevalence - max(platform error, lowest error bar)
            prevalence_range = [
                platform.prevalence_error,
                self.prevalence - max(
                    platform.prevalence_error,
                    max([m['error_bars'][0] for m in mutations]
                ))
            ]

        error = platform.prevalence_error

        # we must allow a little more wiggle for T224i-containing samples
        # because D222S appears quite frequently at low prevalence.
        # TODO: generalize this somehow
        if "T69i" in [str(s) for s in self.sequence.drms]:
            error *= 1.25


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
                'error_bars': [error, error],
                'mutations': mutations, 
                'calls': calls,
                'notes': notes,
                'prevalence_range': prevalence_range
                }
