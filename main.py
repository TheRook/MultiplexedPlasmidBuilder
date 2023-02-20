import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import typer

app = typer.Typer()

@app.command()
def modify_plasmid(
    backbone_file: str = typer.Option(..., "--backbone-file", "-b", help="Path to the GenBank file of the plasmid backbone to modify."),
    target_dir: str = typer.Option(..., "--target-dir", "-t", help="Path to the directory containing the list of target genes to replace. The target genes should be in FASTA format."),
    replacement_dir: str = typer.Option(..., "--replacement-dir", "-r", help="Path to the directory containing the list of replacement genes to use. The replacement genes should be in FASTA format."),
    crispr_file: str = typer.Option(..., "--crispr-file", "-c", help="Path to the GenBank file containing the CRISPR target sequence and PAM."),
    output_file: str = typer.Option(..., "--output-file", "-o", help="Path to the output GenBank file with the modified plasmid.")
):
    """
    Modify a GenBank file of a plasmid backbone by replacing specified target genes with new sequences.

    Parameters:
        backbone_file (str): Path to the GenBank file of the plasmid backbone to modify.
        target_dir (str): Path to the directory containing the list of target genes to replace. The target genes should be in FASTA format.
        replacement_dir (str): Path to the directory containing the list of replacement genes to use. The replacement genes should be in FASTA format.
        crispr_file (str): Path to the GenBank file containing the CRISPR target sequence and PAM.
        output_file (str): Path to the output GenBank file with the modified plasmid.
    """
    # Parse the CRISPR GenBank file to get the target sequence and PAM
    crispr_record = SeqIO.read(crispr_file, "genbank")
    target_sequence = crispr_record.features[0].qualifiers["target"][0]
    pam_sequence = crispr_record.features[0].qualifiers["PAM"][0]

    # Parse the GenBank file of the plasmid backbone
    backbone_record = SeqIO.read(backbone_file, "genbank")
    
    # Loop through each target gene and its replacement
    for target_file, replacement_file in zip(os.listdir(target_dir), os.listdir(replacement_dir)):
        # Parse the target gene and its replacement
        target_record = SeqIO.read(os.path.join(target_dir, target_file), "fasta")
        replacement_record = SeqIO.read(os.path.join(replacement_dir, replacement_file), "fasta")
        
        # Find the location of the target gene in the plasmid backbone
        target_feature = None
        for feature in backbone_record.features:
            if feature.type == "CDS" and feature.qualifiers.get("translation") == target_record.seq.translate():
                target_feature = feature
                break
        
        # Replace the target gene with the replacement sequence
        if target_feature is not None:
            # Find the location of the target sequence in the plasmid backbone
            target_index = backbone_record.seq.find(target_sequence)
            if target_index == -1:
                raise ValueError(f"Could not find target sequence {target_sequence} in plasmid backbone.")
            
            # Find the location of the PAM in the plasmid backbone
            pam_index = backbone_record.seq.find(pam_sequence, target_index+len(target_sequence))
            if pam_index == -1:
                raise ValueError(f"Could not find PAM sequence {pam_sequence} downstream of target sequence {target_sequence} in plasmid backbone.")
            
            # Replace the target sequence with the replacement sequence
            replacement_feature = SeqFeature(
                FeatureLocation(target_index, pam_index+len(pam_sequence)),
                type="CDS",
                qualifiers={
                    "translation": str(replacement_record.seq.translate()),
                    "product": target_feature.qualifiers["product"],
                    "note": target_feature.qualifiers["note"],
                },
            )
            backbone_record.features.remove(target_feature)
            backbone_record.features.append(replacement_feature)
        else:
            raise ValueError(f"Could not find target gene {target_record.id} in plasmid backbone.")
    
    # Write the modified plasmid GenBank file to disk
    with open(output_file, "w") as f:
        SeqIO.write(backbone_record, f, "genbank")

if __name__ == "__main__":
    app()
