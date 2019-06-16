"""
file parser for HPLC data exported as csv
data name and detector identification is based on the MANUALLY GIVEN FILENAME!!  this is crucial since the exported
data doesnt contain any meta-information

the parser is based heavily on the chemstation parser written by Kenneth Nielsen, meant to result in a class structure
that is accessible by the same data treatment code.

author: Anna Winiwarter
date: 14.6.2019
"""

#main feature: Sequence-type class


from __future__ import print_function, unicode_literals, division
import os
import numpy as np


class HplcSequence(object):
    """ The Sequence class for the Chemstation data format

    Parameters:
    as similar to the chemstation data format as possible
        injections (list): List of :class:`~Injection`'s in this sequence
        sequence_dir_path (str): The path of this sequence directory
        metadata (dict): Dict of metadata
    """

    def __init__(self, sequence_dir_path):
            """Instantiate object properties
            Args:
                injections (list): list of Injection objects containing data and metadata for each injection
                sequence_dir_path (str): The path of the sequence
            """
            self.injections = []
            self.sequence_dir_path = sequence_dir_path
            self.metadata = {}
            self._parse()
            # if not self.injections:
            #     msg = 'No injections in sequence: {}'.format(self.sequence_dir_path)
            #     raise NoInjections(msg)
            # self._parse_metadata()

    def _parse(self):
        #needs rewrite because now we dont have a folder for each injection, but only a number of files with the same
        #structures in the filename
        #should fill both the list of injections - need to make a class for that too, and another subclass for
        #the actual data. and that is definitely overkill for the problem at hand, but well...

        sequence_dircontent = os.listdir(self.sequence_dir_path)
        # Put the injection folders in order
        sequence_dircontent.sort()
        previous_filename = "xyz_D"
        for filename in sequence_dircontent:
            # print(filename)
            if filename.endswith(".csv") or filename.endswith(".CSV"):
                if filename[:filename.find("_D")] == previous_filename[:previous_filename.find("_D")]:
                    print("Injection processed already. Continuing.")
                    continue
                else:
                    injection_fullpath = os.path.join(self.sequence_dir_path, filename)
                    self.injections.append(Injection(injection_fullpath))
                    previous_filename = filename


class Injection(object):
    """The Injection class for the Chemstation data format
        kind of the core of this script, collects the raw data from different files with a common filename into one.
        Parameters:
            injection_dirpath (str): The path of the directory of this injection
            metadata (dict): Dict of metadata: sample name, "sample_info" which is also sample name
            detectors (list): list of detectors that are found for this injection
            raw_files (dict): Mapping of detector name (from filename) -> :class:`~HplcFile` objects

        """

    # instantiate the injection object
    def __init__(self, injection_dirpath):
        # print(os.path.dirname(injection_filepath))

        #problem: injection_dirpath is a file not a directory path!!
        print("creating an Injection object.")
        # print(injection_filepath)
        self.injection_dirpath = injection_dirpath
        self.metadata = {}
        self.raw_files = {}
        self.detectors = []
        self._load_raw_spectra(injection_dirpath)

    def _load_raw_spectra(self, injection_filepath):
        """Load all the raw spectra associated with this injection
        creates a dictionary with detector_key as the key (need to read that from the filename!!)
        and a HplcData object (that is similar to the CHFile object from chemstation)
        """
        injection_dirpath = os.path.dirname(injection_filepath)
        filename = os.path.split(injection_filepath)[1]

        for file_ in os.listdir(injection_dirpath):
            if filename in file_:
                print("Creating injection object related to file " + str(filename))
                filepath = os.path.join(injection_dirpath, file_)
                detector = file_[file_.find("_D_")+3:file_.find("_I")] #this here is of course VERY likely to cause errors
                # unless I really, really name files consistently
                self.raw_files[detector] = HplcData(filepath)
                self.detectors.append(detector)
                self.metadata['sample_info'] = "." + file_[:file_.find("_D_")] + "12345"
                self.metadata['sample_name'] = file_[:file_.find("_D_")]

class HplcData(object):
    """
    simplified version of CHFile clas
    this really doesnt need to be a class for what it does, but necessary to result in a file structure as similar to
    chemstation parser as possible
        Attributes:
            values (numpy.array): The intensity values (y-value) or the spectrum. The unit
                for the values is given in `metadata['units']`
            metadata (dict): The extracted metadata #this could be sample name and the detector from the filename or maybe just make attribute .detector?
            filepath (str): The filepath this object was loaded from
        """

    def __init__(self, filepath):
        self.filepath = filepath
        self.values = None #needs to be a numpy array of the y values
        self.times = None #needs to be numpy array of x-values (time)
        self._parse()

    def _parse(self):
        # print(self.filepath)
        file_content = np.genfromtxt(self.filepath, delimiter=",", encoding="UTF16")
        # print(file_content)
        self.content = file_content
        self.values = file_content[:, 1]
        self.times = file_content[:, 0]

