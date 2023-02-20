## Modify Plasmid GenBank File Script
This is a Python script that modifies a GenBank file of a plasmid backbone by replacing specified target genes with new sequences. The script can be run from the command line using the typer package.

### CLI Input
The script accepts the following input parameters:

- --backbone-file: Required. The path to the GenBank file of the plasmid backbone to modify.
- --target-dir: Required. The path to the directory containing the list of target genes to replace. The target genes should be in FASTA format.
- --replacement-dir: Required. The path to the directory containing the list of replacement genes to use. The replacement genes should be in FASTA format.
- --crispr-file: Required. The path to the GenBank file containing the CRISPR target sequence and PAM.
- --output-file: Required. The path to the output GenBank file with the modified plasmid.


### Example Usage
To use the script, navigate to the directory containing the script and run the following command:

`python modify_plasmid.py modify_plasmid --backbone-file /path/to/backbone.gb --target-dir /path/to/targets/ --replacement-dir /path/to/replacements/ --crispr-file /path/to/crispr.gb --output-file /path/to/output.gb`

This command will modify the plasmid GenBank file located at /path/to/backbone.gb by replacing the target genes located in /path/to/targets/ with the replacement genes located in /path/to/replacements/, using the CRISPR target sequence and PAM located in /path/to/crispr.gb, and write the modified plasmid to the file located at /path/to/output.gb.
