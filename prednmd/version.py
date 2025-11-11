"""
Version information for predNMD package
"""

__version__ = "1.0.0"

def get_version_string():
    """Return formatted version string"""
    return f"predNMD v{__version__}"

def get_vcf_annotation_lines(command=None):
    """
    Return VCF header annotation lines for this version
    Following VCF 4.2 specification for meta-information lines
    
    Args:
        command: The command line used to run the software (optional)
    """
    lines = [
        f"##predNMD_version={__version__}"
    ]
    if command:
        lines.append(f"##predNMD_command={command}")
    return lines

def get_table_annotation_lines(command=None):
    """
    Return comment lines for tab-delimited output tables
    
    Args:
        command: The command line used to run the software (optional)
    """
    lines = [
        f"# predNMD version: {__version__}"
    ]
    if command:
        lines.append(f"# Command: {command}")
    return lines
