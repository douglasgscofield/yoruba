yoruba
======

Yoruba is a toolset to query and manipulate BAM files.  Yoruba is under active
development.

Yoruba uses the BamTools C++ API for handling BAM files
(<https://github.com/pezmaster31/bamtools>), and SimpleOpt for handling
command-line options (<http://code.jellycan.com/simpleopt>).

~~~~~~

Usage:   yoruba readgroup [options] <in.bam>
         yoruba kojopodipo [options] <in.bam>

Either command invokes this function.

Add or replace read group information in the BAM file <in.bam>.

Read group information appears in two places in the BAM file:

   (1) the read group dictionary, found in the header, which contains the 
       read group ID and any other information associated with the ID, such 
       as library, sample name, etc., and thus defines the read group;
   (2) the RG tag on each read, which specifies one of the IDs that appear
       in the read group dictionary, and thus declares the read to be part
       of the given read group.

By default, all reads in the BAM file will be given the supplied read group.
If the dictionary already defines a read group with the same ID, its definition
will be replaced with the supplied information.  If the dictionary contains
other read groups, their definitions will remain in the BAM file header but
all reads will be given the supplied read group.

Other behaviour can be specified using --no-replace, --only-replace and 
--clear-read-group.  See table below.

The only argument required to specify a valid read group is --ID or --id.

Options: --ID STR | --id STR                 read group identifier
         --LB STR | --library STR            library
         --SM STR | --sample-name STR        sample name
         --DS STR | --description STR        description
         --DT STR | --date STR               date
         --PG STR | --programs STR           programs used
         --PL STR | --platform STR           sequencing platform
         --PU STR | --platform-unit STR      platform unit
         --PI STR | --predicted-insert STR   predicted median insert size
         --FO STR | --flow-order STR         flow order
         --KS STR | --key-sequence STR       key sequence
         --CN STR | --sequencing-center STR  sequencing center
         --o FILE | -o FILE | --output FILE  output file name [default is stdout]
         --no-replace                        abort if the read group exists
         --only-replace                      replace just this read group
         --clear-read-group                  clear all read group information
         --? | -? | --help                   longer help

         --debug INT     debug info level INT
         --reads INT     process at most this many reads
         --progress INT  print reads processed mod INT

No formatting restrictions are imposed on any of the read group elements. It
is the user's responsibility to ensure that they conform to the SAM definitions
<http://samtools.sourceforge.net/SAM1.pdf> or to any other tool requirements.

If the output file is not specified, then output is written to stdout.

The --no-replace option will abort if the given read group ID is found in the
dictionary, and will only add read group information to reads that don't
don't already have it.

The --only-replace option modifies information for only those reads in the
supplied read group (same ID). Read group information for other reads,
including those without any other read group information, is unchanged.

The --clear-read-group option removes all read group information from all reads.
If specified with options defining a read group, then the read group dictionary
will be cleared prior to defining the new read group.

To summarizing the effects of these options:

                          Read read group (RG) tag status                          
                    ---------------------------------------------                  
                        no RG    |    RG matches   |  RG does not                  
Option                           |       --ID      |  match --ID    RG dictionary  
------------------  ---------------------------------------------  ----------------
                                                                                   
only --ID etc.                 new RG set for all reads             RG added       
                                                                                   
--no-replace         new RG set         abort         no change     RG added; abort
                     from --ID                                      if present     
                                                                                   
--only-replace       no change        no change       no change     RG updated     
                                                                    from options   
                                                                                   
--clear-read-group   no change        RG removed      RG removed    cleared        
  no --ID                                                                          
                                                                                   
--clear-read-group             new RG set for all reads             cleared, then  
  with --ID                                                         RG added       


Kojopodipo is the Yoruba (Nigeria) word for 'group'.

~~~~~~
