yoruba
======

Yoruba is a toolset to query and manipulate BAM files.  Yoruba is under active
development.  Yoruba has an interface similar to [samtools](http://samtools.sourceforge.net)
and some other tools:

    yoruba <command> [options] [<in.bam>] ...

where `<command>` is one of several specific commands.  Thus far, only `readgroup`
is completely implemented.

`readgroup`
: Add or replace read group information

`contents`
: Summarize BAM file contents

Yoruba uses the BamTools C++ API for handling BAM files
(<https://github.com/pezmaster31/bamtools>), and SimpleOpt for handling
command-line options (<http://code.jellycan.com/simpleopt>).


readgroup
---------

    yoruba readgroup [options] [<in.bam>]
    yoruba kojopodipo [options] [<in.bam>]

Add or replace read group information in a BAM file.  *Kojopodipo* is the
Yoruba (Nigeria) verb for 'to group'.  Either command invokes this function.  If
`<in.bam>` is not supplied, input is read from `stdin`.  At most one input BAM
file is allowed.

`yoruba readgroup` is faster and uses less memory than picard `AddOrReplaceReadGroups`.
For a 208GB BAM file containing 10.4M reference sequences and 2.41B reads, 
`AddOrReplaceReadGroups` required ~30 h and ~9 GB RAM to complete, while `yoruba readgroup`
required ~21 h and ~6 GB RAM. 

Read group information appears in two places in a BAM file:

1. the *read group dictionary*, found in the header, which contains definitions
   of individual read groups including the read group ID and any other
   information associated with the ID, such as library, sample name, etc.

2. the *RG tag on each read*, which specifies an ID that appears in the read
   group dictionary, and declares the read to be part of the identified read
   group

By default, all reads in the BAM file will be given the supplied read group.
If the dictionary already defines a read group with the same ID, its definition
will be replaced with the supplied information.  If the dictionary contains
other read groups, their definitions will remain in the BAM file header (if
present) but all reads will be given the supplied read group.

This behaviour can be changed by using the options `--replace` and `--clear`.
See table below.

The only argument required to specify a valid read group is `--ID` or its
synonym `--id`.

| Option                                     | Description |
|--------------------------------------------|-------------|
| `--ID STR` or `--id STR`                   | read group identifier |
| `--LB STR` or `--library STR`              | read group library |
| `--SM STR` or `--sample-name STR`          | read group sample name |
| `--DS STR` or `--description STR`          | read group description |
| `--DT STR` or `--date STR`                 | read group date |
| `--PG STR` or `--programs STR`             | read group programs used |
| `--PL STR` or `--platform STR`             | read group sequencing platform |
| `--PU STR` or `--platform-unit STR`        | read group platform unit |
| `--PI STR` or `--predicted-insert STR`     | read group predicted median insert size |
| `--FO STR` or `--flow-order STR`           | read group flow order |
| `--KS STR` or `--key-sequence STR`         | read group key sequence |
| `--CN STR` or `--sequencing-center STR`    | read group sequencing center |
| `--o FILE` or `-o FILE` or `--output FILE` | output file name [default is stdout] |
| `--replace` STR                            | replace read group STR with --ID
| `--clear`                                  | clear all read group information |
| `--?` or `-?` or `--help`                  | longer help |
| `--debug INT`                              | debug info level `INT` |
| `--reads INT`                              | process at most this many reads |
| `--progress INT`                           | print reads processed mod `INT` |

In the options table, `STR` indicates a string argument, `INT` indicates an
integer value, and `FILE` indicates a filename.

No formatting restrictions are imposed on any of the read group strings. It is
the user's responsibility to ensure that they conform to the SAM definitions
(<http://samtools.sourceforge.net/SAM1.pdf>) or to any other tool requirements.

If the output file is not specified, then output is written to stdout.

The `--replace` option will replace the identified read group to have the name
provided in `--ID`, in both its dictionary entry and on reads.  If only `--ID`
is provided, then the read group is simply renamed.  If any other read group
options are given, then the read group is redefined as well.

The `--clear` option removes all read group information from all reads.  If
specified with options defining a read group, then the read group dictionary
will be cleared prior to defining the new read group.

Only one of these may be supplied at a time.  To summarize the effects of these
options on the read group dictionary and the RG tag on reads:

<table>
<thead>
<tr bgcolor="#e4e4e4">
  <th bgcolor="#e4e4e4">Option</th>
  <th align="center" colspan="3" bgcolor="#e4e4e4">Read Group (RG) tag on reads</th>
  <th bgcolor="#e4e4e4">RG dictionary</th>
</tr>
<tr>
  <th bgcolor="#e4e4e4"></th>
  <th bgcolor="#eeeeee">no RG</th>
  <th bgcolor="#eeeeee">RG matches <code>STR</code></th>
  <th bgcolor="#eeeeee">RG does not match <code>STR</code></th>
  <th bgcolor="#e4e4e4"></th>
</tr>
</thead>
<tbody>
<tr>
  <td>only <code>--ID</code>, etc.</td>
  <td align="center" colspan="3">new RG set for all reads</td>
  <td align="center">RG added</td>
</tr>
<tr>
  <td><code>--replace STR</code></td>
  <td align="center">no change</td>
  <td align="center">RG changed to <code>--ID</code></td>
  <td align="center">no change</td>
  <td align="center">RG <code>STR</code> updated with <code>--ID</code>; replaced if any other RG options</td>
</tr>
<tr>
  <td><code>--clear</code>, no <code>--ID</code></td>
  <td align="center">no change</td>
  <td align="center">RG removed</td>
  <td align="center">RG removed</td>
  <td align="center">cleared</td>
</tr>
<tr>
  <td><code>--clear</code>, with <code>--ID</code></td>
  <td align="center" colspan="3">new RG set for all reads</td>
  <td align="center">cleared, then RG added</td>
</tr>
</tbody>
</table>



contents
--------

    yoruba contents [options] [<in.bam>]
    yoruba inu [options] [<in.bam>]

Summarizes the contents of the BAM file.  *Inu* is the Yoruba (Nigeria) noun
for 'inside'.  Either command invokes this function.  If `<in.bam>` is not
supplied, input is read from `stdin`.  At most one input BAM file is allowed.
No changes to the BAM file are caused by use of this command.

The contents of a BAM file are printed in six sections, the first five comprise
the header and the last is the reads.  The sections in the order described
in the SAM definition (<http://samtools.sourceforge.net/SAM1.pdf>):

1. the *header line* (`@HD`) contains BAM metadata, and containing the tags
   `VN` (required if the header line is present), `SO` and `GO`

2. the *reference sequences* (`@SQ`) describe the reference sequences to which
   the reads in the BAM are alined, and contains the tags `SN` and `LN` (both
   required), and optionally ÀS`, `M5`, `SP`, and `UR`

3. the *read group dictionary* (`@RG`), described under `readgroup` above

4. the *program chain* (`@PG`) describes programs which have manipulated the
   BAM file, and contains the tags ÌD` (required), `PN`, `CL`, `PP` and `VN`.
   
5. *comment lines* (`@CO`) which are individual text lines

6. *reads*, which may be aligned or unaligned.


| Option                     | Description |
|----------------------------|-------------|
| `--refs-to-report` INT     | number of reference sequences to provide details about |
| `--reads-to-report` INT    | number of reads to provide details about [10] |
| `--raw-to-report` INT      | number of characters from the header to print, with `--raw` |
| `--quit`                   | stop after reporting detailed reads, do not continue |
| `--quiet`                  | do not report summaries; still checks header validity and `--raw` |
| `--raw`                    | print raw header information |
| `--?` or `-?` or `--help`  | longer help |
| `--debug INT`              | debug info level `INT` |
| `--reads INT`              | process at most this many reads |
| `--progress INT`           | print reads processed mod `INT` |

In the options table, `INT` indicates an integer value.


