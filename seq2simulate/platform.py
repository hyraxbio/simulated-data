import os

class Platform:

    def __init__(
        self,
        coverage,
        mean_read_length,
        profile=None
    ):
        self.coverage = coverage
        self.profile = profile
        self.mean_read_length = mean_read_length