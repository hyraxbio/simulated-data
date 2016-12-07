import re
from functools import total_ordering

PR, RT, INI = range(3)
locus_names = {
    PR: "PR",
    RT: "RT",
    INI: "INI"
}
loci = {v: k for k,v in locus_names.iteritems()}

drm_regex = re.compile("([A-Z])([0-9]+)([A-Zid])")
wildtype_regex = re.compile("([A-Z])([0-9]+)")
locus_start_positions_tuple = [(INI, 715), (RT, 155), (PR, 56)]
locus_start_positions = {
    d[0]: d[1] for d in locus_start_positions_tuple
}

class MalformattedDrmString(Exception):
    pass

@total_ordering
class Drm:
    """A mutation conferring drug resistance in HIV."""

    def __init__(self, drm_string, locus=None):

        if locus is not None:
            self.locus = locus

        self.insert = False
        self.delete = False
        drm = drm_regex.match(drm_string)
        if drm is None:
            drm = wildtype_regex.match(drm_string)
            if not drm:
                raise MalformattedDrmString("Not a DRM string: " + drm_string)
            # add the wildtype as the mutation if no mutation specified
            drm_string += drm.group(1)
            # check the match
            drm = drm_regex.match(drm_string)
            if not drm:
                raise MalformattedDrmString("Not a DRM string: " + drm_string)

        self.wildtype = drm.group(1)
        if locus is not None:
            self.relative_pos = int(drm.group(2))
            self.absolute_pos = self.relative_pos \
                + locus_start_positions[self.locus]
        else:
            self.absolute_pos = int(drm.group(2))
            found = False
            for region, position in locus_start_positions_tuple:
                if self.absolute_pos > position:
                    self.relative_pos = self.absolute_pos - position
                    found = True
                    break
            if not found:
                raise ValueError('DRM is not in a known gene region')
        self.nucleotide_pos = (self.absolute_pos - 1) * 3    

        self.mutation = drm.group(3)
        if self.mutation == 'i':
            self.insert = True
        elif self.mutation == 'd':
            self.delete = True
        else:
            self.mutation = drm.group(3)

    def locus_str(self):
        return "[" + locus_names[self.locus] + "] " \
            + str(self.relative_pos) \
            + self.mutation

    def __key(self):
        return (self.wildtype, 
                self.absolute_pos, self.mutation)

    def __lt__(self, other):

        return self.__key() < other.__key()

    def __eq__(self, other):

        return self.__key() == other.__key()


    def __hash__(self):
        return hash(self.__key())

    def __str__(self):

        return self.wildtype + str(self.absolute_pos) + self.mutation


    def __repr__(self):
        
        return self.wildtype + str(self.absolute_pos) + self.mutation     






        