import numpy

class MotifProbabilities(object):

    motifs = {
             'GA':   0.125,
             'GG':   0.875,
             }

    def __init__(self):
        for motif, value in self.motifs.iteritems():
            if value == 0.0:
                value = 1e-18
                self.motifs[motif] = value
            setattr(self, motif, value)

    def __str__(self):
        probability_strings = '\n'.join(['{}: {}'.format(i, j) for i, j in self.motifs.iteritems()])
        return '<MotifProbabilities\n{}>'.format(probability_strings)
    
    @property
    def motifs_cumulative(self):
        return self._get_cumulative_p(self.motifs)

    def _get_cumulative_p(self, p):
        """
        Returns cumulatively probabilities from probability dict sorted by
        cumulative probability. This function is useful for sampling mutations.
    
        Args:
            p: sorted cumulative probability dict
        """
     
        motifs = p.keys()
        probabilities = p.values()
        p_cumsum = numpy.array(probabilities).cumsum()
        return dict(zip(p_cumsum, motifs))


class KijakProbabilities(MotifProbabilities):
    """
    This class contains the tetranucleotide hypermutation probabilities from
    Kijak (2008), scaled to reflect a GG:GA ratio of 7.

    Refs:
      Kijak, G. H., Janini, L. M., Tovanabutra, S., Sanders-Buell, E., Arroyo, M.
      A., Robb, M. L., McCutchan, F. E. (2008). Variable contexts and levels of
      hypermutation in HIV-1 proviral genomes recovered from primary peripheral blood
      mononuclear cells. Virology, 376(1), 101-111. doi:10.1016/j.virol.2008.03.017
    """

    motifs = {
             'GA': 0.125100080064,
             'CGGC': 0.0,
             'TGGG': 0.127101681345,
             'AGGG': 0.0850680544436,
             'GGGT': 0.0288230584468,
             'CGGG': 0.0880704563651,
             'GGGA': 0.0601481184948,
             'AGGC': 0.00220176140913,
             'TGGA': 0.0926741393114,
             'CGGA': 0.0681545236189,
             'GGGG': 0.0750600480384,
             'AGGT': 0.043034427542,
             'TGGT': 0.0796637309848,
             'GGGC': 0.00730584467574,
             'CGGT': 0.0504403522818,
             'AGGA': 0.0498398718975,
             'TGGC': 0.0173138510809,
             }

