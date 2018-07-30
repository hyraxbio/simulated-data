import os

import sierra_ws as sierra
import hiv_drms

class Sequence:
    """A sequence from which to simulate reads."""

    def __init__(self, susceptible, resistant):
        self.years_infected = 1
        self.susceptible = susceptible
        self.resistant = resistant
        self.drms = sierra.get_drms(resistant, known_drms=hiv_drms.drms)
        self.pcr_error = False
        self.human_error = False
        self.env_error = False
        self.remove_rt = False
        



