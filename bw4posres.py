#!/usr/bin/env python3.7
"""
================================================================================
 File:        bw4posres.py
 Authors:     Mauricio Esguerra
 Date:        June 23, 2015
 Email:       mauricio.esguerra@gmail.com

 Description:
 With this code we wish to do various task in one module:

 1. Translate pdb to fasta without resorting to import Bio.

 2. Align the translated fasta sequence to a Multiple Sequence Alignment (MSA)
    and place Marks coming from a network of identified conserved
    pair-distances of Venkatakrishnan et al.
    clustalo --profile1=GPCR_inactive_BWtags.aln --profile2=mod1.fasta \
    -o withbwtags.aln --outfmt=clustal --wrap=1000 --force -v -v -v

 3. Translate Marks into properly identified residues in sequence. Notice that
    this depends on a dictionary which uses the Ballesteros-Weinstein numbering.

 4. From sequence ID. pull the atom-numbers of corresponding c-alphas
    in the matched residues.

================================================================================
"""
import os
import subprocess

import settings as s

import logging
logger_bw4posres = logging.getLogger('pymemdyn.bw4posres')

# A larger dictionary of three to one letters can be used.
# For example the dictionary contained at Data/SCOPData.py in a Biopython
# distribution.
# Usually the path where biopython is installed will look something like
# the following path: /lib/python2.7/site-packages/Bio
# For now we are just extending the dictionary to the various possible cases
# for Histidines, that is, HIA, HIC, HID, HIE, HIP, HIQ
protein_letters_3to1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIA': 'H', 'HIC': 'H',
                        'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'HIQ': 'H', 'HIS': 'H',
                        'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
                        'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y',
                        'VAL': 'V'}

bwtags = ['1.46', '1.49', '1.50', '1.53', '1.57', '2.42', '2.43', '2.44', '2.47', '2.50',
          '3.34', '3.36', '3.38', '3.40', '3.44', '3.46', '3.47', '3.51', '4.50', '4.53',
          '4.57', '5.54', '5.57', '5.60', '6.41', '6.44', '6.47', '6.48', '6.51', '7.38',
          '7.39', '7.45', '7.46', '7.47', '7.50', '7.53']

bwpairs = [('1.46', '7.47'), ('1.49', '7.50'), ('1.50', '2.47'), ('1.50', '2.50'),
            ('1.50', '7.46'), ('1.53', '2.47'), ('1.57', '2.44'), ('2.42', '3.46'),
            ('2.43', '7.53'), ('2.50', '7.46'), ('3.34', '4.53'), ('3.34', '4.57'),
            ('3.36', '6.48'), ('3.38', '4.50'), ('3.38', '4.53'), ('3.40', '6.44'),
            ('3.44', '5.54'), ('3.47', '5.57'), ('3.51', '5.57'), ('3.51', '5.60'),
            ('5.54', '6.41'), ('6.47', '7.45'), ('6.51', '7.38'), ('6.51', '7.39')]

class Run(object):
    """
    A pdb file is given as input to convert into one letter sequence
    and then align to curated multiple sequence alignment and then
    assign Ballesteros-Weinstein numbering to special positions.
    """
    def __init__(self, pdb, **kwargs):
        """
        The init method is a kind of constructor, called when an instance
        of the class is created. The method serves to initialize what you
        want to do with the object.
        """
        self.pdb = pdb
        self.own_dir = kwargs.get("own_dir") or ""
        self.clustal_bin = s.CLUSTAL_BIN
        self.repo_dir = kwargs.get("repo_dir") or s.ROOT_DIR + "/templates"
        self.loggerbw = logging.getLogger('pymemdyn.bw4posres.Run')

    def pdb2fas(self):
        """
        From pdb file convert to fasta sequence format without the use of
        dependencies such as BioPython. This pdb to fasta translator
        checks for the existance of c-alpha residues and it is
        based on their 3-letter sequence id.
        """
        fastaseq = open(os.path.join(self.pdb.split(".")[0] + ".fasta"), "w")
        self.loggerbw.debug('writing now: '+os.path.join(self.pdb.split(".")[0] + ".fasta"))
        fastaseq.write(">")
        fastaseq.write("{0}\n".format(self.pdb))
        result = []
        for line in open(self.pdb, "r"):
            atoms = [a for a in line if line[0:6] == "ATOM  " and line[13:16] == "CA "]
            seqnam = ''.join(atoms[17:20])
            if seqnam != '':
                seq = protein_letters_3to1[seqnam]
                result.append(seq)

        lines = []
        numcol = 70
        resultasstr = ''.join(result)
        for i in range(0, len(resultasstr), numcol):
            lines.append(resultasstr[i:i+numcol])

        fastaseq.write("{0}".format("\n".join(lines)))
        fastaseq.write("\n")
        fastaseq.close()

    def clustalalign(self):
        """
        Align the produced fasta sequence with clustalw to assign
        Ballesteros-Weinstein marks.
        """
        profile1 =  self.repo_dir + "/GPCR_inactive_BWtags.aln"
        profile2 = os.path.join(self.pdb.split(".")[0] + ".fasta")
        bwtagged = os.path.join(self.pdb.split(".")[0] + "_bw" + ".aln")

        command = [
            self.clustal_bin,
            "-profile1=" + profile1,
            "-profile2=" + profile2,
            "-outfile=" + bwtagged,
            "-output=clustal"]

        proc = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        proc.communicate()

        return True

    def getcalphas(self):
        """
        Pulls out the atom numbers of c-alpha atoms. Restraints are
        placed on c-alpha atoms.
        """
        calphaspdb = open(os.path.join(self.pdb.split(".")[0] + "_CA.pdb"), "w")
        for line in open(self.pdb, "r"):
            atoms = [a for a in line if line[0:6] == "ATOM  " and line[13:16] == "CA "]
            calphaspdb.write("{0}".format("".join(atoms)))

    def makedisre(self):
        """
        Creates a disre.itp file with atom-pair id's to be restrained
        using and NMR-style Heaviside function based on
        Ballesteros-Weinstein tagging.

        the disre.itp contains the reference distances of pairs of C-alpha at conserved GPCR positions in inactive structures
        used for Ballesteros-Weinstein restrains at the final stage of equilibration
        BW only works for monomers (utils.py line:207)
        """
        disre = open('disre.itp', 'w')
        readdisre =  open(self.repo_dir + "/disre.itp", 'r')
        disre.write("; This file provides distance restraints based on an online example\n")
        disre.write("; by David Van der Spoel for alpha-helices.\n")
        disre.write("; Notice that for now we are using a fixed value for up2 of 1.2.\n")
        disre.write("; low is the average pair distance for the set of 19 inactive structures\n")
        disre.write("; of Ramakrishnan et al. minus one sd, up1 is the average plus one sd.\n")
        disre.write("; ai aj type index type' low up1 up2 fac\n")
        disre.write("[ distance_restraints ]\n")

        # This part pulls out the two last lines from the clustalw alignment
        # and puts each in a list.
        # list1 has the Ballesteros-Weinstein positions for the restraints.
        # list2 has the final alignment of the sequence, corresponding to the
        # pdb that was given.
        list1 = []
        list2 = []
        bwtagged = os.path.join(self.pdb.split(".")[0] + "_bw" + ".aln")
        for line in open(bwtagged, "r"):
            fields = line.split()
            try:
                if fields[0] == 'bw':
                    list1.append(fields[1])
                if fields[0] == self.pdb:
                    list2.append(fields[1])
            except IndexError:
                continue

        list1s = ''.join(list1)
        list1l = list(list1s)
        list2s = ''.join(list2)
        list2l = list(list2s)

        # The two lists are combined into a list of tuples using the zip builtin
        # and then the dashes (sequence gaps) are parsed out of the
        # corresponding well aligned sequence to be able to count without gaps.
        # Finally the residue numbers corresponding to the A mask are returned.
        ziplist = list(zip(list1l, list2l))
        nodash = []
        for index, string in enumerate(ziplist):
            if ziplist[index][1].isalpha():
                nodash.append(ziplist[index])

        resid = []
        for index, string in enumerate(nodash):
            if nodash[index][0] == 'A':
                resid.append(index+1)

        # To map residue id's to c-alpha atom numbers I open the parsed
        # pdb which contains only c-alphas and look for the existence of
        # the located residue numbers in the resid list.
        caid = []
        onlyca = open(os.path.join(self.pdb.split(".")[0] + "_CA.pdb"), "r")
        firstline = onlyca.readline()
        offset =  int(int(firstline[22:26]) - 1)
        renumbered = [x + offset for x in resid]
        logger_bw4posres.debug('renumbered: {}'.format(renumbered))
        # Take care of the possible offset of sequence numbers not starting
        # at 1.
        for line in onlyca:
            if int(''.join(line[22:26])) in renumbered:
                caid.append(line[7:11])

        # These two lists must have the same dimension:
        logger_bw4posres.debug("The following two lists must have the same dimension:")
        logger_bw4posres.debug('{} \t {}'.format(len(bwtags), len(caid)))
        
        # The mapping (rosetta-stone) between BW id numbers and c-alpha atom
        # numbers is made in this list.
        bw2calpha = list(zip(bwtags,caid))

        bwatom1 = []
        bwatom2 = []
        for index, value in enumerate(bwpairs):
            bwatom1.append(value[0])
            bwatom2.append(value[1])

        caatom1 = []
        for i in range(0, 24):
            for j in range (0,36):
        # TODO: make sure that a pdb with gaps will also go through.
        # right now a pdb with gaps on the sequential numbering of
        # residues will have a wrong mapping and not even an error message
        # is issued.
        # Issue reported by user Laurens Kooijman
                if bwatom1[i] == bw2calpha[j][0]:
                    caatom1.append(bw2calpha[j][1])

        caatom2 = []
        for i in range(0, 24):
            for j in range (0,36):
                if bwatom2[i] == bw2calpha[j][0]:
                    caatom2.append(bw2calpha[j][1])

        for i in range(0,24):
            logger_bw4posres.debug('{}, {}, {}'.format(caatom1[i], caatom2[i], bwpairs[i]))

        # This part finally writes to the disre.itp file all the information it needs
        # for the atom pairs after they have been mapped from Ballesteros-Weinstein
        # into c-alpha atom numbers in the original pdb given as input.
        for index, string in enumerate(readdisre):
            fields = string.split()
            fields.insert(0, caatom2[index])
            fields.insert(0, caatom1[index])
            #It would be better to have first 9 fields be float
            #instead of str.
            fields.insert(9, '\n')
            disre.write('{0}'.format('   '.join(fields)))
        disre.close()

        return True