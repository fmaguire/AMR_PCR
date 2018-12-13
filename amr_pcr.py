#!/usr/bin/env python

import pandas as pd
import json
import re


def read_primers(primer_tsv_fp):
    """
    Parse the primer data into a clean pandas table
    """
    # tidy up and extract data using tabula
    primers = pd.read_csv(primer_tsv_fp)
    primers['Sequence_tidy'] = primers['Sequence'].str.replace(' ', '').str.replace("3", "").str.replace("5", "").str.replace("-", "").str.replace("’", "").str.replace("'", "")

    # removing MRSA as not sure what this is
    primers = primers[primers['Antimicrobial'] != 'MRSA']

    # tidy up filled columns
    primers['CARD_name'] = primers['Gene/Codon no.']

    na_group = primers[primers['CARD_name'].isna()].index
    primers.loc[na_group, 'CARD_name'] = primers.loc[na_group, 'Group']

    na_group = primers[primers['CARD_name'] == 'all'].index
    primers.loc[na_group, 'CARD_name'] = primers.loc[na_group, 'Group']

    primers['CARD_name'] = primers['CARD_name'].str.replace(' group', '')


    # equivalent names to CARD naming scheme
    equivalent_names = {'catA1': 'catI',
                        'vat': 'vatA',
                        'VanA': 'vanA',
                        'VanB': 'vanB',
                        'aac(3)-IV': 'AAC(3)-IV',
                        'ant(2")-I': "ANT(2'')-Ia",
                        'aac(3)-II': 'AAC(3)-II',
                        'aph(3’)-III': "APH(3')-III",
                        'aph(3’)-II': "APH(3')-II",
                        'aph(3’)-I': "APH(3')-I",
                       "aac(6') Ib-cr": "AAC(6')-Ib-cr",
                       "strA":"APH(3'')-Ib",
                       "strB":"APH(6)-Id",
                       'aadE': "ANT(6)-I",
                       "tetA": "tet(A)",
                       "tetB": "tet(B)",
                       "tetC": "tet(C)",
                       "tetD": "tet(D)",
                        "tetE": "tet(E)",
                        "tetG": "tet(G)",
                        "Tet (H)": "tet(H)",
                        'tet(W)': "tetW",
                        'tet(M)': 'tetM',
                        'tet(O)': 'tetO',
                        'tet(S)': 'tetS',
                        'Tet (T)': 'tetT',
                        'Tet (Z)': 'tet(Z)',
                        'Tet (W)': 'tetW',
                        'Tet (31)': 'tet(31)',
                        'tet(32)': 'tet32',
                        'Tet (34)': 'tet34',
                        'Van X': 'vanX',
                        'CTX-M9': 'CTX-M-9',
                        'CTX-M1': 'CTX-M-1',
                        'CTX-M2': 'CTX-M-2',
                        'M-All': 'CTX-M',
                        'qnrA': 'QnrA',
                        'qnrB': 'QnrB',
                        'qnrC': 'QnrC',
                        'qnrD': 'QnrD',
                        'qnrS': 'QnrS',
                        'qepA': 'QepA',
                        'gyrA E. Coli': 'gyrA Escherichia',
                        'parC E. Coli': 'parC Escherichia',
                        "ANT(2'')-I": "ANT(2'')-Ia"}


    primers['CARD_name'] = primers['CARD_name'].replace(equivalent_names)

    return primers


class CARD():
    """
    Class to contain CARD metadata
    """
    def __init__(self, card_fp):
        with open(card_fp) as fh:
            self.card = json.load(fh)

        print("CARD version " + self.card['_version'])
        del self.card['_version']
        del self.card['_comment']
        del self.card['_timestamp']


        # so we get all the CARD genes that are meant to be detected by a given
        # primer set we need to set up some parsing rules

        # ACC-1 might just be ACC, DHA, VEB

        # primer should match any CARD entry with a dash suffix e.g. TEM-2
        self.dash_indices = ['TEM',
                             'FOX',
                             'SHV',
                             'IMP',
                             'SPM',
                             'VIM',
                             'KPC',
                             'NDM',
                             'CTX-M']

        # primer should match the name followed by a number e.g. qnrA4
        self.num_indices = ['cmlA',
                            'QnrA',
                            'QnrB',
                            'QnrD',
                            'QnrS',
                            'QepA',
                            'aadA2']


        # get all with following letters AACX-III{a,b,c,...}
        self.letter_indices = [
                   "AAC(3)-II",
                   "APH(3')-III",
                   "APH(3')-II",
                   "APH(3')-I",
                   "ANT(6)-I",
                   "vanX"]

        # treating tet as specific hits
        #specific single genes for exact matches
        self.specific = ['ACC-1', 'DHA-1', 'OXA-48',
                         'QnrC',"AAC(3)-IV", "ANT(2'')-Ia",
                         'catI', 'floR', "AAC(6')-Ib-cr",
                         "APH(3'')-Ib", "APH(6)-Id",
                         'sul1', 'sul2', 'sul3', 'tet(A)', 'tet(B)',
                         'tet(C)', 'tet(D)', 'tet(E)', 'tet(G)',
                         'tet(H)', 'tet(K)', 'tet(L)',
                         'tetM', 'tetO', 'tetS', 'tetT', 'tetW',
                         'tet(Z)', 'tet(31)', 'tet32', 'tet(33)', 'tet34',
                         'tet(39)', 'vanA', 'vanB',
                         'ErmA', 'ErmB', 'ErmC', 'ErmF', 'vatB', 'vatD', 'vatE',
                         'vgaA', 'vgaB', 'vgbA', 'vgbB', 'vatA', 'ErmE',
                         'CMY-1',
                         'CMY-2', 'VEB-1',
                         'CTX-M-1',
                         'CTX-M-2',
                         'CTX-M-9']

        # name must contain species and name
        self.species_specific = ['gyrA Salmonella',
                    'parC Salmonella',
                    'gyrA Escherichia',
                    'parC Escherichia']


        # Unsure if CMY and CTX are meant to be more than singles
        # Unknown what CTX-M-ALL is meant to hit but listed as every
        # BIC doesn't seem to be in CARD by that name

    def get_aro_for_name(self, name):
        """
        Return the set of AROs that are meant to hit for a given name
        """
        aro_names = []
        aros_accessions = []

        if name == 'BIC':
            return [], []

        for key, entry in self.card.items():

            # if its one of the exact match AROs
            if name in self.specific:
                if entry['ARO_name'] == name:
                    aro_names.append(entry['ARO_name'])
                    aros_accessions.append(entry['ARO_accession'])

            # if its one of the species specific matches that should have
            # a couple of hits
            elif name in self.species_specific:
                matches = []
                for part in name.split():
                    for entry_part in entry['ARO_name'].split():
                        if part == entry_part:
                            matches.append(True)
                if all(matches) and len(matches) == len(name.split()):
                    aro_names.append(entry['ARO_name'])
                    aros_accessions.append(entry['ARO_accession'])

            elif name in self.num_indices:
                m = re.search(re.escape(name) + '\d+$', entry['ARO_name'])
                if m is not None:
                    aro_names.append(entry['ARO_name'])
                    aros_accessions.append(entry['ARO_accession'])

            elif name in self.letter_indices:
                m = re.search(re.escape(name) + '[a-z,A-Z]+$', entry['ARO_name'])
                # if the string ends in letters m will be a match
                if m is not None:
                    aro_names.append(entry['ARO_name'])
                    aros_accessions.append(entry['ARO_accession'])

            elif name in self.dash_indices:
                m = re.search(re.escape(name) + '-\d+$', entry['ARO_name'])
                # if the string ends in digits m will be a Match object, or None otherwise.
                if m is not None:
                    aro_names.append(entry['ARO_name'])
                    aros_accessions.append(entry['ARO_accession'])

        if len(aro_names) == 0:
            print("ARO match missing ", name)
            assert False
        return aro_names, aros_accessions
