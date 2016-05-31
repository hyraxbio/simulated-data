import os

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