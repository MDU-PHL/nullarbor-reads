# `nullarbor-reads`

## Goal
The goal of this script is to produce a `nullarbor` input file. It is
designed to be quite flexible.

## Running

### Creating input file by regex pattern

Assuming one has a folder with read files, which may or may not be organised
into subfolders according to isolate/sample, and you want all the reads that
match a regular expression 'myreads[0-9]{4}'. So, each read file starts with
`myreads` and is followed by exactly four numbers, you would do the following:

        nullarbor-reads --seq_path /path/to/seqs/folder --id_pattern myreads[0-9]{4} input.tab

The `input.tab` is the only argument for `nullarbor-reads`, and is the output
filename, where the `nullarbor` input file will be saved to.

`--seq_path` only needs to be defined if outside of the sequence folder. The
default for `--seq_path` is `'.'`.

If you want to see what is going on, run it with `--verbose`:

        nullarbor-reads --id_pattern myreads[0-9]{4} --verbose input.tab

### Creating input file by idfile

Assuming you have a tab-delimited file (TSV) with one ore more columns, and one
column has the ID of the isolates of interest, one can run

        nullarbor-reads --seq_path /path/to/seqs/folder --idfile isolates.txt input.tab

If the ID column is not the first (but the 3rd, for instance), and there is a header row,
just use the following:

        nullarbor-reads --seq_path /path/to/seqs/folder --idfile isolates.txt --col_number 3 --header_true input.tab

### If you want to search deeper sub-folders levels

Assuming your sequence folder has more than one level of subfolders:

        folder/
            subfolder1/
                sub-subfolder1/
                    seq_reads1.fastq.gz
            subfolder2/

You can increase the maximum level to search with the `--level` flag:

        nullarbor-reads --seq_path /path/to/seqs/folder --idfile isolates.txt --level 2 input.tab

### If you want to exclude files/folders with certain keywords

Assuming you have multiple sequence files for each isolate/sample, but you want
to exclude some:

        nullarbor-reads --seq_path /path/to/seqs/folder --idfile isolates.txt --exclude old --exclude CLIPPED input.tab

In the above, any files/subfolders with that have `old` or `CLIPPED` in the name
will be excluded. You can add as many `--exclude` as you need.

### If your read files don't have a traditional extension, or there are a mix of extensions

By default, `nullarbor-reads` will search for any files with the following four
extensions: fastq, fq, fastq.gz, fq.gz. If you want to add to this list just use
the following:

        nullarbor-reads --seq_path /path/to/seqs/folder --idfile isolates.txt --alt_extension fa input.tab

### If your PE reads are not named with 'R1' and 'R2'

By default `nullarbor-reads` expects that read files will be annotated with 'R1'
and 'R2' to distinguish among the pair of files for a single sample. Similarly to
`--exclude` you can add as many new strings to distinguish among read files in
a pair as you want. So, if you had read pairs that were separated by `read1` and
`read2`, you would do the following:

        nullarbor-reads --seq_path /path/to/seqs/folder --idfile isolates.txt --read1_pat read1 --read2_pat read2 input.tab

## TODO

    1. Add some logic to resolve conflict. If there is more than a single pair of files
    which one should be chosen.
        --- One possibility is to use last date modified.
        
