"""
utils.py — Biophysical Utilities

This module contains shared utility functions used across the HelixZero-CMS engine. 
These functions evaluate structural properties of the siRNA sequences that influence
RISC loading, off-target thermodynamic stability, and innate immunogenicity.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)


def calculate_gc_percentage(sequence: str) -> float:
    """
    Calculates the GC content of an RNA sequence.
    
    The thermodynamic stability of the siRNA duplex is heavily influenced by 
    the ratio of Guanine-Cytosine pairs (which have 3 hydrogen bonds) to 
    Adenine-Uracil pairs (which have 2). Optimal GC content (typically 30-55%) 
    ensures that the duplex is stable enough in serum but loose enough to be 
    unwound by the RNA-induced silencing complex (RISC).
    
    Args:
        sequence (str): The nucleotide string (RNA/DNA) to analyze.
        
    Returns:
        float: The GC percentage (0.0 to 100.0). Returns 0.0 if sequence is empty.
    """
    if not sequence:
        logger.debug("Empty sequence provided for GC calculation.")
        return 0.0
        
    upper_seq = sequence.upper()
    gc_count = upper_seq.count("G") + upper_seq.count("C")
    return (gc_count / len(sequence)) * 100.0


def has_internal_palindrome(sequence: str, half_length: int = 4) -> bool:
    """
    Detects internal palindromic regions within an RNA sequence.
    
    Internal palindromes (like a hairpin loop) cause the single-stranded siRNA 
    to fold back on itself. If the Antisense strand forms a strong internal 
    secondary structure, it will sterically hinder binding to the target mRNA 
    and block the Ago2 slicing action. This function scans for such formations.
    
    Args:
        sequence (str): The RNA sequence to analyze.
        half_length (int): The required length of the palindromic stem. Defaults to 4.
        
    Returns:
        bool: True if a reverse-complement palindrome of at least `half_length` exists, False otherwise.
    """
    if len(sequence) < half_length * 2:
        return False
        
    # Translation table to generate the reverse complement for RNA
    trans_table = str.maketrans("AUGC", "UACG")
    
    for i in range(len(sequence) - 2 * half_length + 1):
        # Extract the half-segment, reverse it, and complement it
        reverse_complement = sequence[i : i + half_length][::-1].translate(trans_table)
        
        # Check if the reverse complement appears anywhere downstream of the current segment
        if reverse_complement in sequence[i + half_length :]:
            logger.debug(f"Palindromic region found in sequence: {reverse_complement}")
            return True
            
    return False
