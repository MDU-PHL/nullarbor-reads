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
            self.reads['R1'] = os.path.join( self.work_path, read_path )
        else:
            self.reads['R1'] = os.path.join( self.work_path, read_path )

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

    def is_complete(self):
        '''This function returns True if there are two read sets associated with
        this SequenceID object, and False otherwise.
        This will help when printing out to the input file, and identifying any
        sequence IDs that are not complete.
        This could happen if the user has a typo in their file or naming of the
        sequence IDs, or if the sequence reads do not exist.
        '''
        if( length( self.reads.keys() ) == 2):
            return True
        else:
            return False

class SequenceCollection:
    '''
    A class that holds a collection of SequencesID objects.
    It has the capacity to create SequenceID objects, produce summary information
    on the objects, and print out the Nullarbor input file.
    '''
    def __init__(self):
        ''' Initiate the object with a dictionary to hold SequenceID objects
        '''
        self.sequences = {}

    def add_from_file(self, fn, sequence_path, col_number = 0, header = False, alternate_extension = None):
        '''
        If the user has a file with a list of Sequence IDs that they wish to
        include in their input file, this is the easiest way to generate a
        Nullarbor input file.

        Behaviour:
            1. Reads the list of IDs from the file --- the file should be
            formatted as a single column text file, with one ID per column ---
            if more than one column provided, the user should provide a column
            number to identify the column ID. The user should also specify if
            there is a column header or not (default no column header)
            2. For each successfuly read ID, initiate a SequenceID object.
            3. For each SequenceID object, search the sequence path for
            compatible read sets. Here, we have to allow for some leeway in
            searching --- possible files could have have the following extensions
            fastq, fq, fastq.gz, fq.gz. However, we do allow for an alternate
            extension provided by the user.
            4. If the collection of reads found is different from two, try to
            figure out which ones to keep. There is no simple logic here. We
            could decide to keep the largest file, or the latest file, or one
            that has a certain extension. First, the program will try to group
            the files into pairs that have similar styles, and then apply some
            selection criteria.
            5. Once two read sets have been identified, then place them in the
            SequenceID object
            '''

    def add_from_folder(self, sequence_path, alternate_extension = None):
        '''
        Similar to add_from_file.

        Behaviour:
            1. Search sequence_path for unique pairs of files (possibly in
            unique subfolders)
            2. Find ID from file pairs, and then group pairs by ID. Make
            SequenceID objects
            3. If more than a pair of reads per ID, then apply rationale
            in add_from_file to figure out which pair to keep.
            4. Add suitable pair of reads to SequenceID object
        '''

    def make_nullarbor_input(self, outfn):
        '''
        Using outfn as a filename, print Nullarbor input.
        Behaviour:
            1. Check all SequenceID objects to ensure they have the
            appropriate number of reads. Create list of sequence IDs
            to include in file.
            2. Print summary information to screen to warn user of any
            missing information
            3. Open outfn, generate a single string with all input lines, write
            to outfn, and close outfn.
            4. Issue some warning to the screen to tell the user of what has
            happened.
        '''

class SearchReads():
    '''
    This class will be responsible for searching and matching reads into
    suitable pairs.
    '''
    def __init__(self, alt_extension = None):
        '''
        When initialising this object, we only need to set up a dictionary to
        keep all the read sets, grouped by ID, and some pattern to find read
        files.
        '''
        self.read_groups = {}
        self.extension_pat = ("fastq", "fq", "fastq.gz", "fq.gz")
        if (alt_extension != None):
            self.extension_pat = self.extension_pat + (alt_extension,)

    def find_all_reads(self, sequence_path, level = 1):
        '''
        This function finds all reads in the sequence path going down into
        X levels in the sub-directory path.
        '''
        reads_path = os.path.abspath(sequence_path)
        read_files = []
        num_seps = reads_path.count( os.path.sep ) # will need this below to
        # figure out where I am in the levels
        for dirname, subdirs, files in os.walk( reads_path ):
            # here I am using the number of path separators to identify
            # at which level I am in the directory structure. So, the
            # root directory will have num_seps, and each additional level
            # will add an additional separator. If the number of additional
            # separators, beyond that in the sequence_path are found, then
            # we skip that folder.
            if( ( dirname.count( os.path.sep ) - num_seps ) > level ):
                continue
            for f in files:
                print(f)
                if f.endswith( self.extension_pat ):
                    read_files.append( os.path.join( dirname, f ) )
        return read_files

    def search_by_id(self, id_list, sequence_path):
        '''
        Find sequences that are matched by ID. This will be used by
        add_by_id function above.
        The output is a dictionary, kept in self.read_groups, which has
        sequenceID as a key, and different lists of paired reads as the values.
        '''

    def search_by_path(self, sequence_path, id_pattern):
        '''
        Find sequences, and then infer ID. This will be used by
        add_by_path function above.
        Here, the user has to give a regular expression to search for and
        extract the ID from.
        '''

@click.command()
@click.option("--seq_path", default = ".", help = "Path to sequence reads (default = '.')")
@click.option("--idfile", default = None, help = "File with sequence IDs")
@click.option("--col_number", default = 0, help = "If idfile has more than one column, which column has the sequence IDs")
@click.option("--header_true", is_flag=True, help = "If idfile has a header row")
@click.option("--alt_extension", default = None, help = "If sequence files don't have one of the following extensions: fastq, fq, fastq.gz, fq.gz, specify alternate extension. Search will add this extension to the default ones, so mixed extensions is allowed (e.g., 'txt' ).")
@click.option("--id_pattern", default = None, help = "If searching by folder rather than using sequence ID file, what is the regex pattern to identify the ID from a file name.")
def nullarbor_reads(seq_path, idfile, col_number, alt_extension, id_pattern, header_true):
    '''
    This is the main function
    '''
    new_search = SearchReads()
    tmp = new_search.find_all_reads(seq_path)
    print(tmp)
    return

if __name__=="__main__":
    nullarbor_reads()
