# By Scott Teresi
# The below code contains my class for genes and transposons
from collections import defaultdict
#-------------------------------------
class Genic_Element(object):
    def __init__(self, number, chromosome, start, stop ):
        self.number = number
        self.chromosome = chromosome
        self.start = int(start)
        self.stop = int(stop)

    def getNumber(self):
        return self.number

    def getChromosome(self):
        return self.chromosome

    def getStart(self):
        return self.start

    def getStop(self):
        return self.stop

    def getLength(self):
        return self.length

class TE(Genic_Element):
    def __init__(self, number, chromosome, start, stop, length, te_type, family):
        super().__init__(number, chromosome, start, stop)
        self.te_type = te_type
        self.family = family
        self.length = length

    def getTe_Type(self):
        return self.te_type

    def getFamily(self):
        return self.family


class Gene(Genic_Element):
    def __init__(self, number, chromosome, start, stop, maker_name):
        super().__init__(number, chromosome, start, stop)
        self.maker_name = maker_name

        # I am going to create a left, center, and right density for each TE element, including classes and families
        # This is where things get obnoxiously complicated.
        # Classes are my two distinctions (referred to in my code as type and fam)
        # Classes:
        self.DNA_left = 0
        self.DNA_intra= 0
        self.DNA_right = 0

        self.LTR_left = 0
        self.LTR_intra= 0
        self.LTR_right = 0

        self.Unknown_left = 0
        self.Unknown_intra= 0
        self.Unknown_right = 0

        self.LINE_left = 0
        self.LINE_intra= 0
        self.LINE_right = 0
        #-----------------
        # families

        self.MULE_left = 0
        self.MULE_intra= 0
        self.MULE_right = 0

        self.Gypsy_left = 0
        self.Gypsy_intra= 0
        self.Gypsy_right = 0

        self.Unknown_fam_left = 0
        self.Unknown_fam_intra= 0
        self.Unknown_fam_right= 0

        self.CMC_EnSpm_left = 0
        self.CMC_EnSpm_intra= 0
        self.CMC_EnSpm_right = 0

        self.Copia_left = 0
        self.Copia_intra= 0
        self.Copia_right = 0

        self.None_left = 0
        self.None_intra= 0
        self.None_right = 0

        self.hAT_left = 0
        self.hAT_intra= 0
        self.hAT_right = 0

        self.PIF_Harbinger_left = 0
        self.PIF_Harbinger_intra= 0
        self.PIF_Harbinger_right = 0

        #------------------------
        # Proximities
        self.prox_left = None
        self.prox_right = None
        #------------------------
        # TE_Count
        self.they_are_inside = 0


    def getMaker_Name(self):
        return self.maker_name
    def getFrag_Name(self):
        return self.Frag_Name
    def getleft(self):
        return self.left_density
    def getright(self):
        return self.right_density
    def getintra(self):
        return self.intra_density
    def getProxRight(self):
        return self.prox_right
    def getProxLeft(self):
        return self.prox_left


if __name__ == '__main__':
    pass
