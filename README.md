## Modify Plasmid GenBank File Script
This project simplifies the process of building multiplexed-plasmids and packaging them up into a single annoated genbank file that can be easily used on benchling or another workbench.  Although this task is commonly perfomred by hand in a browser or desktop application, performing senstive tasks like this with python allows us to be eaxct with our changes and to verify the changes.

### CLI Input
The script accepts the following input parameters:

- --backbone-file: Required. The path to the GenBank file of the plasmid backbone to build from.
- --target-dir: Required. The path to the directory containing the list of target genes to replace.
- --replacement-dir: Required. The path to the directory containing the list of replacement genes to use.
- --crispr-file: Required. The path to the GenBank file containing the desired multiplexed-CRISPR variant.
- --output-file: Required. The path to the final plasmid.


### Example Usage
To use the script, navigate to the directory containing the script and run the following command:

`python modify_plasmid.py modify_plasmid --backbone-file /path/to/backbone.gb --target-dir /path/to/targets/ --replacement-dir /path/to/replacements/ --crispr-file /path/to/crispr.gb --output-file /path/to/output.gb`

This command will modify the plasmid GenBank file located at /path/to/backbone.gb by replacing the target genes located in /path/to/targets/ with the replacement genes located in /path/to/replacements/, using the CRISPR target sequence and PAM located in /path/to/crispr.gb, and write the modified plasmid to the file located at /path/to/output.gb.
