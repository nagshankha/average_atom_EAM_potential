#!/usr/bin/env python
# -*- coding: utf-8 -*-
#################
from inspect import currentframe, getframeinfo
import numpy as np
import time
import os
#################

"""
Module containing class MixPotentials

Author: Created by Shankha Nag
"""
## -------------------------------------------------------------------------- ##


class MixPotentials():

    """ 
class MixPotentials: Creating an equivalent potential from potentials of different elements by averaging over concentrations

Embedding potential: F_A = sum_i(conc[i] * F[i])
Electron density:    rho_A = sum_i(conc[i] * rho[i])
Pair potential:      phi_r_Aj = sum_i(conc[i] * phi_r_ij),     (sum_i = summation over i, conc = respective concentration)
                     phi_r_AA = sum_j(conc[j] * phi_r_Aj)

The counter i runs over all the elements that are considered for mixing potential.

    """

    ## ------------------------------------------------------------------ ##
    ## Constructors/Destructors                                           ##
    ## ------------------------------------------------------------------ ##

    # The Atomtype_conc must be of type 'OrderedDict' from 'collections' class
    def __init__(self, pot_file, Atomtype_conc, new_elname):
        """ Instantiation:

Member variables initialized during Instantiation:

PotFileName (type: string) :: Path to the potential file
Elname_conc (type: dict) :: Stores the element names for averaging as keys and the corresponding concentrations as values
new_elname (type: string) :: Name of the new element that needs to be created

The setInputs() method is used for initializing the members Elname_conc and new_elname to check the type compatibility. 


 """

        ## Members ---------------------- ##

        # string PotFileName
        self.PotFileName = pot_file

        # dict Elname_conc
        # string new_elname

        self.setInputs(Atomtype_conc, new_elname)

        pass

    def __del__(self):
        """ pass """
        pass

    ## ------------------------------------------------------------------ ##
    ## Methods                                                            ##
    ## ------------------------------------------------------------------ ##

    # public:

    def readPotential(self):
        """ Reads the Potential file and gathers the information regarding embedding potential, electron densities and pair potentials """

        if (not isinstance(self.PotFileName, str)):
            raise TypeError('class MixPotentials (line ' + str(getframeinfo(
                currentframe()).lineno) + '): Potential Filename must be a string')

        elif (len(self.PotFileName) <= 10):
            raise RuntimeError('class MixPotentials (line ' + str(getframeinfo(currentframe()).lineno) +
                               '): The Potential filename should have more than 10 characters (for .eam.alloy extension at the end)')
        elif (self.PotFileName[-10:] != '.eam.alloy'):
            raise RuntimeError('class MixPotentials (line ' + str(getframeinfo(
                currentframe()).lineno) + '): Potential filename should have a .eam.alloy extension')

        if (not os.path.exists(self.PotFileName)):
            raise RuntimeError('class MixPotentials (line ' + str(getframeinfo(
                currentframe()).lineno) + '): Potential file does not exist in the specified path')

        f = open(self.PotFileName, 'r')

        initial_string = 'Mixed Interatomic potential developed by Shankha Nag where New Element \'' + \
            self.new_elname + '\' is made of {'

        for i in list(self.Elname_conc.keys()):
            initial_string += '\'' + i + '\': ' + \
                str(self.Elname_conc[i]) + ', '
        initial_string = initial_string[:-2] + '}\n'

        initial_string += ('Based Interatomic Potential file from: ' + self.PotFileName + '\n' +
                           'Mixing algorithm implemented via a python code on ' + time.strftime("%d/%m/%Y") + ' at ' + time.strftime("%H:%M:%S") + '\n')

        for i in range(0, 3):
            next(f)

        line = next(f)
        nEl_basePot = int(line.split()[0])

        if (len(set(line.split()[1:]).intersection(set(self.Elname_conc.keys()))) != len(set(self.Elname_conc.keys()))):
            raise RuntimeError('class MixPotentials (line ' + str(getframeinfo(currentframe(
            )).lineno) + '): The provided potential file does not have the required atom types')
        else:
            ElIndex = [line.split()[1:].index(x)
                       for x in list(self.Elname_conc.keys())]
            initial_string += str(len(line.split()[1:]) + 1) + ' ' + ' '.join(
                line.split()[1:]) + ' ' + self.new_elname + '\n'

        line = next(f)
        Nrho = int(line.split()[0])
        Nr = int(line.split()[2])
        initial_string += line

        mass = []
        F = []  # Embedding Function
        rho = []  # Density Function
        FLineCount = 0
        rhoLineCount = 0
        mass_F_rho_string = ''

        for i in range(0, nEl_basePot):
            F.append([])
            rho.append([])

            for j in range(0, 3):

                if (j == 0):
                    line = next(f)
                    mass.extend([float(line.split()[1])])
                    if (i == 0):
                        if (len(line.split()[1].split('.')[1].split('e')[0]) != len(line.split()[1].split('.')[1])):
                            formNotMass = 'e'
                            massPrecision = len(
                                line.split()[1].split('.')[1].split('e')[0])
                        elif (len(line.split()[1].split('.')[1].split('E')[0]) != len(line.split()[1].split('.')[1])):
                            formNotMass = 'E'
                            massPrecision = len(
                                line.split()[1].split('.')[1].split('E')[0])
                        else:
                            formNotMass = 'f'
                            massPrecision = len(line.split()[1].split('.')[1])
                    mass_F_rho_string += line

                elif (j == 1):
                    dataCount = 0
                    while (dataCount < Nrho):
                        if (i == 0):
                            FLineCount += 1
                        line = next(f)
                        F[i].extend([float(k) for k in line.split()])
                        dataCount += len(line.split())
                        mass_F_rho_string += line
                    if (i == 0):
                        if (len(line.split()[0].split('.')[1].split('e')[0]) != len(line.split()[0].split('.')[1])):
                            formNotF = 'e'
                            FPrecision = len(line.split()[0].split('.')[
                                             1].split('e')[0])
                        elif (len(line.split()[0].split('.')[1].split('E')[0]) != len(line.split()[0].split('.')[1])):
                            formNotF = 'E'
                            FPrecision = len(line.split()[0].split('.')[
                                             1].split('E')[0])
                        else:
                            formNotF = 'f'
                            FPrecision = len(line.split()[0].split('.')[1])

                elif (j == 2):
                    dataCount = 0
                    while (dataCount < Nr):
                        if (i == 0):
                            rhoLineCount += 1
                        line = next(f)
                        rho[i].extend([float(k) for k in line.split()])
                        dataCount += len(line.split())
                        mass_F_rho_string += line
                    if (i == 0):
                        if (len(line.split()[0].split('.')[1].split('e')[0]) != len(line.split()[0].split('.')[1])):
                            formNotRho = 'e'
                            rhoPrecision = len(
                                line.split()[0].split('.')[1].split('e')[0])
                        elif (len(line.split()[0].split('.')[1].split('E')[0]) != len(line.split()[0].split('.')[1])):
                            formNotRho = 'E'
                            rhoPrecision = len(
                                line.split()[0].split('.')[1].split('E')[0])
                        else:
                            formNotRho = 'f'
                            rhoPrecision = len(line.split()[0].split('.')[1])

        phi_r = []  # Pair Interaction times r
        phi_rLineCount = 0
        phi_r_string = ''

        for i in range(0, int(nEl_basePot * (nEl_basePot + 1) / 2)):
            phi_r.append([])

            dataCount = 0
            while (dataCount < Nr):
                if (i == 0):
                    phi_rLineCount += 1
                line = next(f)
                phi_r[i].extend([float(k) for k in line.split()])
                dataCount += len(line.split())
                phi_r_string += line
            if (i == 0):
                if (len(line.split()[0].split('.')[1].split('e')[0]) != len(line.split()[0].split('.')[1])):
                    formNotPhi_r = 'e'
                    phi_rPrecision = len(
                        line.split()[0].split('.')[1].split('e')[0])
                elif (len(line.split()[0].split('.')[1].split('E')[0]) != len(line.split()[0].split('.')[1])):
                    formNotPhi_r = 'E'
                    phi_rPrecision = len(
                        line.split()[0].split('.')[1].split('E')[0])
                else:
                    formNotPhi_r = 'f'
                    phi_rPrecision = len(line.split()[0].split('.')[1])

        f.close()

        formNot = [formNotMass, formNotF, formNotRho, formNotPhi_r]

        precision = [massPrecision, FPrecision, rhoPrecision, phi_rPrecision]

        lineCount = [FLineCount, rhoLineCount, phi_rLineCount]

        return (nEl_basePot, Nrho, Nr, mass, F, rho, phi_r, ElIndex, lineCount, formNot, precision, initial_string, mass_F_rho_string, phi_r_string)

    #######################################################################################################################################

    def compute(self, outdir=None, s_add=None):
        """ Calculates the embedding potential, electron densities and pair potentials for the newly created average element """

        nEl_basePot, Nrho, Nr, mass, F, rho, phi_r, ElIndex, lineCount, formNot, precision, initial_string, mass_F_rho_string, phi_r_string = self.readPotential()

        conc = [self.Elname_conc[i] for i in list(self.Elname_conc.keys())]

        mass_F_rho_A_string = ''

        mass_F_rho_A_string += '100 ' + '{0:.{1}{2}}'.format(np.sum(np.array(conc) * np.array(
            [mass[x] for x in ElIndex])), precision[0], formNot[0]) + ' 100 fcc\n'

        F_A = np.zeros(Nrho)
        rho_A = np.zeros(Nr)

        for i in range(0, len(ElIndex)):

            F_A += conc[i] * np.array(F[ElIndex[i]])
            rho_A += conc[i] * np.array(rho[ElIndex[i]])

        for i in np.array_split(F_A, lineCount[0]):

            mass_F_rho_A_string += ' '.join(
                ['%.{0}{1}'.format(precision[1], formNot[1]) % x for x in i]) + '\n'

        for i in np.array_split(rho_A, lineCount[1]):

            mass_F_rho_A_string += ' '.join(
                ['%.{0}{1}'.format(precision[2], formNot[2]) % x for x in i]) + '\n'

        phi_r_A_string = ''
        phi_r_AR = []

        for i in range(0, nEl_basePot):
            phi_r_AR.append(np.zeros(Nr))
            for j in range(0, len(ElIndex)):

                if (i >= ElIndex[j]):
                    phi_r_AR[i] += conc[j] * \
                        np.array(phi_r[int(i * (i + 1) / 2) + ElIndex[j]])
                else:
                    phi_r_AR[i] += conc[j] * \
                        np.array(
                            phi_r[int(ElIndex[j] * (ElIndex[j] + 1) / 2) + i])

        phi_r_AA = np.sum(np.repeat(np.array(conc).reshape(len(conc), 1), len(
            phi_r_AR[0]), axis=1) * np.array(phi_r_AR)[ElIndex], axis=0)

        for i in range(0, nEl_basePot):
            for j in np.array_split(phi_r_AR[i], lineCount[2]):

                phi_r_A_string += ' '.join(
                    ['%.{0}{1}'.format(precision[3], formNot[3]) % x for x in j]) + '\n'

        for i in np.array_split(phi_r_AA, lineCount[2]):

            phi_r_A_string += ' '.join(
                ['%.{0}{1}'.format(precision[3], formNot[3]) % x for x in i]) + '\n'

        string_2_write = initial_string + mass_F_rho_string + \
            mass_F_rho_A_string + phi_r_string + phi_r_A_string

        self.getResults(string_2_write, outdir, s_add)

    #######################################################################################################################################

    def getResults(self, string_2_write, outdir, s_add=None):
        """ Writes the new potential file and updates the PotFileName member with the filename of the newly created potential file """

        s = self.PotFileName
        if s_add is None:
            s_add = 'MixedPotential_' + \
                '-'.join(list(self.Elname_conc.keys()) +
                         [self.new_elname]) + '.eam.alloy'
        elif not isinstance(s, str):
            raise ValueError('"s_add" must be a string')

        if outdir == None:
            newPotFile = '/'.join(s.split('/')[:-1]) + [s_add if (
                len(s.split('/')) == 1) else ('/' + s_add)][0]
        else:
            if outdir[-2:] == '//':
                raise RuntimeError('Directory path cannot end with "//"')
            elif outdir[-1] != '/':
                outdir += '/'

            if not os.path.isdir(outdir):
                os.makedirs(outdir)

            newPotFile = outdir + s_add

        self.PotFileName = newPotFile

        fw = open(newPotFile, 'w')

        fw.write(string_2_write)

        fw.close()

    ########################################################################################################################################

    def setInputs(self, Atomtype_conc, new_elname):
        """ Initializes the members Elname_conc and new_elname after checking the type compatibility of the values provided by the user """

        if (isinstance(Atomtype_conc, dict)):
            self.Elname_conc = Atomtype_conc
        else:
            raise TypeError('class MixPotentials (line ' + str(getframeinfo(currentframe()).lineno) +
                            '): Atom type name and corresponding concentration must be provided in dictionary')

        for i in list(self.Elname_conc.keys()):
            if (not isinstance(self.Elname_conc[i], float)):
                raise TypeError('class MixPotentials (line ' + str(getframeinfo(
                    currentframe()).lineno) + '): Concentrations of each atom type must be in fraction')

        if not np.isclose(sum([self.Elname_conc[i] for i in list(self.Elname_conc.keys())]), 1):
            raise RuntimeError('class MixPotentials (line ' + str(getframeinfo(
                currentframe()).lineno) + '): Sum of concentrations of all atom types must be unity')

        if (isinstance(new_elname, str)):
            self.new_elname = new_elname
        else:
            raise TypeError('class MixPotentials (line ' + str(getframeinfo(
                currentframe()).lineno) + '): Name of New Homogeneous Element must be a string')


## -------------------------------------------------------------------------- ##
            



if __name__ == '__main__':
    test = MixPotentials()





