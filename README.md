yoruba
======

Yoruba is a toolset to query and manipulate BAM files.  Yoruba is under active
development.  Yoruba has an interface similar to [samtools](http://samtools.sourceforge.net)
and some other tools:

    yoruba <command> [options] [<in.bam>] ...

where `<command>` is one of several specific commands.  Thus far, only `readgroup`
is completely implemented.

| `yoruba` command   | Action |
|--------------------|--------|
| `yoruba readgroup` | Add or replace read group information in a BAM file |

Yoruba uses the BamTools C++ API for handling BAM files
(<https://github.com/pezmaster31/bamtools>), and SimpleOpt for handling
command-line options (<http://code.jellycan.com/simpleopt>).

readgroup or kojopodipo
-----------------------

    yoruba readgroup [options] [<in.bam>]
    yoruba kojopodipo [options] [<in.bam>]

Add or replace read group information in a BAM file.  *Kojopodipo* is the
Yoruba (Nigeria) word for 'group'.  Either command invokes this function.  If
`<in.bam>` is not supplied, input is read from `stdin`.  At most one input BAM
file is allowed.

Read group information appears in two places in a BAM file:

1. the *read group dictionary*, found in the header, which contains definitions
   of individual read groups including the read group ID and any other
   information associated with the ID, such as library, sample name, etc.

2. the *RG tag on each read*, which specifies an ID that appears in the read
   group dictionary, and declares the read to be part of the identified read
   group

By default, `readgroup` will give all reads the supplied read group.  If the
read group dictionary already defines a read group with the same ID, its
definition will be removed and replaced with a definition containing the read
group information given in the command-line options.  If the dictionary defines
other read groups, these definitions will remain in the BAM file header, but
all reads in the BAM file will still be given the new read group.

This behaviour can be changed by using the options `--no-replace`, `--only-replace` 
and `--clear`.  See table below.

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
| `--no-replace`                             | abort if the read group exists |
| `--only-replace`                           | replace just this read group |
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

The `--no-replace` option will abort if the given read group ID is found in the
dictionary, and will only add read group information to reads that don't
don't already have it.

The `--only-replace` option modifies information for only those reads in the
supplied read group (same ID). Read group information for other reads,
including those without any other read group information, is unchanged.

The `--clear` option removes all read group information from all reads.
If specified with options defining a read group, then the read group dictionary
will be cleared prior to defining the new read group.

Only one of these options may be supplied at a time.  To summarize the effects
of these options on the read group dictionary and the RG tag on reads:

<table>
<thead>
<tr bgcolor="#e4e4e4">
  <th bgcolor="#e4e4e4">Option</th>
  <th align="center" colspan="3" bgcolor="#e4e4e4">Read Group (RG) tag on reads</th>
  <th bgcolor="#e4e4e4">RG dictionary</th>
</tr>
<tr>
  <th bgcolor="#e4e4e4"></th>
  <th bgcolor="#eaeaea">no RG</th>
  <th bgcolor="#eaeaea">RG matches <code>--ID</code></th>
  <th bgcolor="#eaeaea">RG does not match <code>--ID</code></th>
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
  <td><code>--no-replace</code></td>
  <td align="center">new RG set</td>
  <td align="center">abort</td>
  <td align="center">no change</td>
  <td align="center">RG added; abort if present</td>
</tr>
<tr>
  <td><code>--only-replace</code></td>
  <td align="center">no change</td>
  <td align="center">no change</td>
  <td align="center">no change</td>
  <td align="center">RG updated from options</td>
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


### Informal performance comparison

| BAM file                           | `yoruba readgroup` | picard `AddOrReplaceReadGroups` |
|------------------------------------|--------------------|---------------------------------|
| 208G with 2.41B reads and 10.4M references | ~23:00:00, 4GB RAM | 31:23:43, >9GB RAM |


