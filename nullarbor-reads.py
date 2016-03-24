#!/usr/bin/env python2
'''
A script to generate a Nullarbor input file from a random folder of
sequence reads
'''

import click
import os
import re

WINDOW_SIZE = 50

class SequenceReads:
    '''A class to store the sequence IDs'''
    def __init__(self, id ):
        '''
        The class will be initiated by creating a dictionary to store reads
        paths to and the ID information
        '''
        self.id =  id
        self.reads = {}

    def add_read(self, read_path, read_type ):
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

        self.reads[ read_type ] = read_path

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
        if( len( self.reads.keys() ) == 2):
            return True
        else:
            return False

class SequenceCollection:
    '''
    A class that holds a collection of SequencesID objects.
    It has the capacity to create SequenceID objects, produce summary information
    on the objects, and print out the Nullarbor input file.
    '''
    def __init__(self, verbose, is_SE):
        ''' Initiate the object with a dictionary to hold SequenceID objects
        '''
        self.sequences = {}
        self.verbose = verbose
        self.is_SE = is_SE

    def add_from_searchreads(self, searchreads, resolve_conflict_by, read1_pat, read2_pat):
        ''' This function will take the output from SearchReads object, and organise into
        read pairs, assuming paired end reads by default. As Nullarbor might accept single
        end reads in the future, I am coding to include this as a possibility.

        Behaviour:
        1. IF PE and less than 2 read files, issue warning, skip isolate, and go to end.
        2. If PE, check that there are at least two files, and that there is at
        least one that matches one of elements in tuple read1_pat and at least one that
        matches the read2_pat elements.
        3. If PE and only 2 read files, and each one matches one of the read1_pat and the other
        matches the other read2_pat, then assign, and skip to end.
        4. Else If PE and no read file matches either one or both of the read file name patterns,
        issue a warning, and skip.
        5. Else If PE and more than one read file matches either one or both of the read file patterns,
        then cluster reads from different read names by string distance.
            -- If after clustering there is only one complete pair, accept that pair, and move to end.
            --- Else if after clustering there are more than one complete pair, then apply
                conflict_resolution_by to pick most likely pair.
            --- Else IF after clustering there are no pairs, then issue warning, skip isolate, and go to end.
        '''
        try:
            read_groups = searchreads.read_groups
        except RuntimeError:
            print("When trying to add_from_searchreads expected object to be of class SearchReads.")
        isolates = read_groups.keys()
        if( self.is_SE ):
            pass
        else:
            for i in isolates:
                reads = read_groups[i]
                n_reads = len( reads )
                if( n_reads == 2):
                    #import pdb; pdb.set_trace()
                    read1 = []
                    read2 = []
                    for r in reads:
                        if( any( pat in r for pat in read1_pat) ):
                            read1.append(r)
                        elif( any( pat in r for pat in read2_pat) ):
                            read2.append(r)
                        else:
                            continue
                    if( len( read1 ) == 1 and len( read2 ) == 1 ):
                        self.sequences[ i ] = SequenceReads( id = i )
                        self.sequences[ i ].add_read(read1[ 0 ], 'R1')
                        self.sequences[ i ].add_read(read2[ 0 ], 'R2')
                    else:
                        if( self.verbose ):
                            print( "One or more files did not match any of read number patterns for isolate {}!".format( i ))
                        continue
                if( n_reads < 2) :
                    if( self.verbose ):
                        print( "Isolate {} HAD LESS THAN 2 reads, and is BEING SKIPPED!".format( i ))
                    continue



    def make_nullarbor_input(self, outfn ):
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
        try:
            out_obj = open( outfn, 'w' )
        except IOError:
            print( "Was not able to open {} to write output.".format( outfn ))
        if( self.verbose ):
            print( '#' * WINDOW_SIZE )
            print( 'Starting to create nullarbor output file...')
        outstring = ''
        successfully_written = 0
        for i in self.sequences:
            tmp = self.sequences[i]
            if ( tmp.is_complete() ):
                tmp = tmp.print_nullarbor()
                outstring = outstring + tmp
                successfully_written = successfully_written + 1
                if( self.verbose ):
                    print( tmp.strip() )
        out_obj.write( outstring )
        out_obj.close()
        if( self.verbose ):
            print( "Found {} isolates, and successfully wrote {} to file {}.".format( len( self.sequences ), successfully_written, outfn ))
        print( "Successfully wrote Nullabor input file. Happy Nullarboring!" )



class SearchReads():
    '''
    This class is responsible for searching and matching reads to IDs.
    '''
    def __init__(self, exclude = None, alt_extension = None, verbose = False):
        '''
        When initialising this object, we only need to set up a dictionary to
        keep all the read sets, grouped by ID, and some pattern to find read
        files.
        '''
        self.read_groups = {}
        self.extension_pat = ("fastq", "fq", "fastq.gz", "fq.gz")
        if (alt_extension != None):
            self.extension_pat = self.extension_pat + (alt_extension,)
        self.exclude = exclude
        self.verbose = verbose

    def find_all_reads(self, sequence_path, level = 1):
        '''
        This function finds all reads in the sequence path going down into
        X levels in the sub-directory path.
        '''
        if( self.verbose ):
            print( '#' * WINDOW_SIZE )
            print("Searching for reads...")
        reads_path = os.path.abspath(sequence_path)
        read_files = []
        num_seps = reads_path.count( os.path.sep ) # will need this below to
        # figure out where I am in the levels
        for dirname, subdirs, files in os.walk( reads_path ):
            # here I am using the number of path separators to identify
            # at which level I am in the directory structure. So, the
            # root directory will have num_seps, and each additional level
            # will add an additional separator. If the number of additional
            # separators is beyond that found in the sequence_path are found, then
            # we skip that folder.
            if( ( dirname.count( os.path.sep ) - num_seps ) > level ):
                continue
            if( any( ex in dirname for ex in self.exclude) ):
                if( self.verbose ):
                    print( "Directory {} was EXCLUDED because it matched one of the excluded patterns".format( dirname ))
                continue
            for f in files:
                if (self.exclude !=None and any( ex in f for ex in self.exclude) ):
                    if( self.verbose ):
                        print( "File {} was EXCLUDED because it matched one of the excluded patterns".format( f ))
                    continue
                if f.endswith( self.extension_pat ):
                    read_files.append( os.path.join( dirname, f ) )
        if( self.verbose ):
            print( "Found {} read files in path {}.".format( len( read_files ), reads_path ))
            print( "" )

        return read_files

    def search_by_id(self, id_list, sequence_path, level, header_true, col_number):
        '''
        Find sequences that are matched by ID. This will be used by
        add_by_id function above.
        The output is a dictionary, kept in self.read_groups, which has
        sequenceID as a key, and different lists of paired reads as the values.

        Behaviour:
            1. Search sequence_path for read files (possibly in
            unique subfolders). Apply any exclude statements to remove
            files and subdirectories from searching
            2. Find ID from file using regex expression provided by the user
            3. Orgainse read files in a dictionary with ID as key
        '''
        sequence_files = self.find_all_reads( sequence_path, level )
        if( self.verbose ):
            print( '#' * WINDOW_SIZE )
            print( "Organizing reads found according to IDS in file {}".format( id_list ))
        # load isolate IDs
        id_index = col_number - 1 # assume that user gives column numbers in 1-base format
        try:
            fn_obj = open( id_list )
        except IOError:
            print("Could not open {}".format( id_list ))
        file_sep = '\t'
        total_ids = 0
        if ( header_true ):
            # skip header if one is present
            fn_obj.next()
        for l in fn_obj:
            tmp = l.strip().split( file_sep )
            tmp = tmp[ id_index ]
            self.read_groups[ tmp ] = []
            total_ids = total_ids + 1
        if( self.verbose ):
            print( "Successfully opened file {}, and found {} IDs".format( id_list, total_ids ))
            print( "" )
        fn_obj.close()
        #import pdb; pdb.set_trace()
        # match reads to keys
        id_keys = self.read_groups.keys()
        not_found = True
        for f in sequence_files:
            for k in id_keys:
                if( re.search( k, f ) ):
                    self.read_groups[k].append( f )
                    not_found = False
                    break
                else:
                    continue
            if( not_found and self.verbose ):
                print( "File {} did not match any id".format( f ))
                not_found = True

    def search_by_path(self, sequence_path, id_pattern, level):
        '''
        Find sequences, and then infer ID. This will be used by
        add_by_path function above.
        Here, the user has to give a regular expression to search for and
        extract the ID from.
        Behaviour:
            1. Search sequence_path for read files (possibly in
            unique subfolders). Apply any exclude statements to remove
            files and subdirectories from searching
            2. Read the list of IDs from the file --- the file should be
            formatted as a single column text file, with one ID column ---
            if more than one column provided, the user should provide a column
            number to identify the column ID. The user should also specify if
            there is a column header or not (default no column header)
            3. For each successfuly read ID, add to read_groups dictionary as key.
            4. Cross reference found read files to keys, add matching read files
            to the list of reads associated with an ID.
        '''
        sequence_files = self.find_all_reads( sequence_path, level )
        if( self.verbose ):
            print( '#' * WINDOW_SIZE )
            print( "Organizing reads found according to pattern {}".format( id_pattern ))
        # searching through files to figure out which match our search criteria
        # and grouping the files by isolate ID, if the ID can be found
        for f in sequence_files:
            try:
                seq_id = re.search( id_pattern, f ).group()
                try:
                    self.read_groups[seq_id].append( f )
                except:
                    self.read_groups[seq_id] = []
                    self.read_groups[seq_id].append( f )
            except:
                if (self.verbose ):
                    print( "File {} DID NOT INCLUDE the specified pattern {}".format( f, id_pattern ))
                continue
        if( self.verbose ):
            print( "" )

    def summary(self):
        '''
        Generating some summary information about the read sets that were found
        '''
        min_no_reads = 2
        max_no_reads = 2
        if( self.verbose ):
            print( '#' * WINDOW_SIZE )
            print( 'Generating some summary infomation about found reads...' )
            print( "Found {} IDs".format( len( self.read_groups.keys() ) ) )
        for isolate in self.read_groups:
            tmp = len( self.read_groups[isolate] )
            if( tmp != 2 and self.verbose ):
                print( "The isolate {} appears to have less/more read sets than expected: {}".format( isolate, self.read_groups[isolate] ) )
                print( "Perhaps use --exclude XXX to reduce the number of reads?")
            if( tmp > max_no_reads ):
                max_no_reads = tmp
                continue
            elif( tmp < min_no_reads ):
                min_no_reads = tmp
                continue
            else:
                continue
        print( "Maximum number of reads for an id was {}".format( max_no_reads ) )
        print( "Minimum number of reads for an id was {}".format( min_no_reads ) )


@click.command()
@click.option("--seq_path", default = ".", help = "Path to sequence reads (default = '.')")
@click.option("--idfile", default = None, help = "File with sequence IDs")
@click.option("--col_number", default = 1, help = "If idfile has more than one column, which column has the sequence IDs")
@click.option("--header_true", is_flag = True, help = "If idfile has a header row")
@click.option("--alt_extension", default = None, help = "If sequence files don't have one of the following extensions: fastq, fq, fastq.gz, fq.gz, specify alternate extension. Search will add this extension to the default ones, so mixed extensions is allowed (e.g., 'txt' ).")
@click.option("--id_pattern", default = None, help = "If searching by folder rather than using sequence ID file, what is the regex pattern to identify the ID from a file name.")
@click.option("--level", default = 1, help = "How many subdirectory levels to go down searching for read files (default = '1')")
@click.option("--exclude", multiple = True, default = None, help = "Any text that can be used to exclude unwanted reads (e.g., from a particular subdirectory, or with a particular pattern. Can be used multiple times (e.g., --exclude subdir1 --exclude old))")
@click.option("--verbose", is_flag = True, help = "Be verbose.")
@click.option("--is_SE", is_flag = True, help = "For future compatibility with Nullarbor accepting single-end reads (default PE reads)")
@click.argument("OUTFILE")
def nullarbor_reads(seq_path, idfile, col_number, alt_extension, id_pattern, header_true, level, exclude, verbose, is_se, outfile):
    '''
    This is the main function
    '''
    read1_pat = ('R1',)
    read2_pat = ('R2',)
    new_search = SearchReads( exclude, alt_extension, verbose )
    if( idfile != None ):
        if( verbose ):
            print( '#' * WINDOW_SIZE )
            print( "Searching {} with ids found in file {}... ".format( seq_path, idfile ) )
            print( '' )
        new_search.search_by_id( idfile, seq_path, level, header_true, col_number)
    elif( id_pattern != None ):
        if( verbose ):
            print( '#' * WINDOW_SIZE )
            print( "Searching {} with pattern {}... ".format( seq_path, id_pattern ) )
            print( '' )
        new_search.search_by_path( seq_path, id_pattern, level )
    else:
        raise RuntimeError("Neither idfile nor id_pattern were specified. One of these needs to be specified.")
    new_search.summary()
    new_collection = SequenceCollection( verbose, is_se )
    new_collection.add_from_searchreads( new_search, resolve_conflict_by = None, read1_pat = read1_pat, read2_pat = read2_pat )
    new_collection.make_nullarbor_input( outfile )
    return

if __name__=="__main__":
    nullarbor_reads()
