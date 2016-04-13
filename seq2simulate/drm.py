import re
from functools import total_ordering

PR, RT, INI = range(3)
gene_names = {
    PR: "PR",
    RT: "RT",
    INI: "INI"
}
drm_regex = re.compile("([A-Z])([0-9]+)([A-Zid])")
gene_start_positions = {PR: 56, RT: 155, INI: 715}

class MalformattedDrmString(Exception):
    pass

@total_ordering
class Drm:

    def __init__(self, drm_string, gene):

        self.gene = gene
        self.insert = False
        self.delete = False
        drm = drm_regex.match(drm_string)
        if drm is None:
            raise MalformattedDrmString("Not a DRM string: " + drm_string)

        self.wildtype = drm.group(1)
        self.relative_pos = int(drm.group(2))
        self.nucleotide_pos = (
                self.relative_pos + gene_start_positions[self.gene] - 1
            ) * 3
        self.mutation = drm.group(3)
        if self.mutation == 'i':
            self.insert = True
        elif self.mutation == 'd':
            self.delete = True
        else:
            self.mutation = drm.group(3)

    def __key(self):
        return (self.gene, self.wildtype, self.relative_pos, self.mutation)

    def __lt__(self, other):

        return self.__key() < other.__key()

    def __eq__(self, other):

        return self.__key() == other.__key()


    def __hash__(self):
        return hash(self.__key())

    def __str__(self):

        return gene_names[self.gene] + ":" + \
        self.wildtype + str(self.relative_pos) + self.mutation


    def __repr__(self):
        
        return gene_names[self.gene] + ":" + \
        self.wildtype + str(self.relative_pos) + self.mutation     






        