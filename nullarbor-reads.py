'''
A script to generate a Nullarbor input file from a random folder of
sequence reads
'''

import click
import os
import re

class SequenceID:
    '''A class to store the sequence IDs'''
    def __init__(self, id, work_path):
        '''
        The class will be initiated by creating a dictionary to store reads
        paths to and the ID information
        '''
        self.id =  id
        self.work_path = work_path # keep track of the work path to ensure
                                    # outfile has the full path
        self.reads = {}

    def add_read(self, read_path):
        '''
        This function adds a read path to the collection
        Behaviour should be:
            Get a read filename found in a folder that matched the same
            ID as the sequence ID.
            Check if it is R1 or R2, create an appropriate entry in the
            dictionary
        Assumptions:
            Read sets always have R1/R2 (but the software allows for r1/r2)
            If there are more than one R1/R2 pair in a folder matching the
            correct ID, then it must be filtered someplace else. There can
            only be two read sets.
        '''
        if (re.search( 'R1', read.path.upper() )):
            self.work_path['R1'] = os.path.join( self.work_path, read_path )
        else:
            self.work_path['R1'] = os.path.join( self.work_path, read_path )

    def print_nullarbor(self, sep = '\t', end_line = '\n'):
        '''
        This function outputs a string in Nullarbor input format. The string
        has the ID, R1, and R2 all separated by tabs, and then a new line
        character at the end.
                sep = '\t' in case separator changes at some point
                end_line = '\n' allow for different end of line characters at some point
        '''
        nullarbor_string = self.id + sep + self.reads['R1'] + sep + self.reads['R2'] + end_line
        return( nullarbor_string )
