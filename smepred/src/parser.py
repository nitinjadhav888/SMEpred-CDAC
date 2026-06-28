"""
parser.py — Sequence Input Parser

Handles the ingestion and normalization of mRNA/gene inputs.
Ensures that all downstream ML models and combinatorial engines 
operate on a sanitized, strictly-RNA string format, preventing 
nucleotide mismatches or length violations during feature extraction.
"""

import re
import logging
from pathlib import Path
from typing import Union

logger = logging.getLogger(__name__)


def _normalize_nucleotides(raw_sequence: str) -> str:
    """
    Strips non-alphabetic characters, enforces uppercase, and converts DNA 'T' to RNA 'U'.
    
    Why: Biological datasets often contain line breaks, whitespace, or are provided as 
    DNA sequences instead of RNA transcripts. Normalization is strictly required to 
    prevent structural biophysics calculations from crashing when encountering 'T' 
    or unexpected whitespace tokens.
    
    Args:
        raw_sequence (str): The raw input sequence.
        
    Returns:
        str: The normalized RNA sequence.
        
    Raises:
        ValueError: If illegal characters remain after normalization.
    """
    # Remove any non-alphabetic characters (including whitespace/newlines) and uppercase
    clean_seq = re.sub(r"[^A-Za-z]", "", raw_sequence).upper()
    
    # Strictly enforce RNA format
    clean_seq = clean_seq.replace("T", "U")
    
    valid_nucleotides = set("AUGC")
    invalid_chars = set(clean_seq) - valid_nucleotides
    
    if invalid_chars:
        logger.error(f"Invalid nucleotides detected: {invalid_chars}")
        raise ValueError(
            f"Input sequence contains unexpected characters: {invalid_chars}. "
            "Only standard A, U, G, C (or T for DNA input) nucleotides are allowed."
        )
        
    return clean_seq


def _extract_first_fasta_sequence(fasta_text: str) -> str:
    """
    Extracts the first sequence from a multiline FASTA formatted string.
    
    Why: Users frequently copy/paste raw FASTA blocks directly from NCBI or Ensembl. 
    This parser isolates the raw string array from the metadata header ('>') block.
    
    Args:
        fasta_text (str): The raw FASTA formatted text.
        
    Returns:
        str: The raw, un-normalized string data of the first sequence.
        
    Raises:
        ValueError: If no sequence data is found below the FASTA header.
    """
    lines = fasta_text.strip().splitlines()
    sequence_lines = []
    is_reading_sequence = False
    
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            # If we were already reading a sequence, hitting a second '>' means we stop.
            if is_reading_sequence:
                break
            is_reading_sequence = True
        elif is_reading_sequence:
            sequence_lines.append(line)
            
    if not sequence_lines:
        logger.error("FASTA extraction failed: No sequence lines found below header.")
        raise ValueError("No sequence data found in the provided FASTA input.")
        
    return "".join(sequence_lines)


def load_sequence(source: Union[str, Path]) -> str:
    """
    Loads and normalizes an mRNA or gene sequence from a file path or raw string.
    
    Acts as the entry point for all target sequences. Automatically detects whether 
    the input is a file path, raw text, or a FASTA block, and routes to the 
    appropriate extractor.
    
    Args:
        source (Union[str, Path]): A file path (.fa, .fasta, .txt) or an inline string.
        
    Returns:
        str: A validated, uppercase RNA sequence containing only A, U, G, C.
        
    Raises:
        ValueError: If the file path is unreadable, FASTA parsing fails, illegal characters are found, or the sequence is too short.
    """
    try:
        path = Path(str(source))
        if path.suffix.lower() in (".fa", ".fasta", ".fna", ".txt") and path.exists():
            logger.info(f"Loading sequence from file path: {path}")
            raw_text = path.read_text(encoding="utf-8")
        else:
            raw_text = str(source)
    except Exception as e:
        logger.error(f"Failed to read source input: {e}")
        # Default to treating it as a string if Path coercion completely fails on weird inputs
        raw_text = str(source)

    if not raw_text.strip():
        logger.error("Empty input provided.")
        raise ValueError("Sequence input cannot be empty.")

    # Detect FASTA format vs plain raw sequence
    if raw_text.lstrip().startswith(">"):
        raw_nucleotides = _extract_first_fasta_sequence(raw_text)
    else:
        raw_nucleotides = raw_text

    normalized_sequence = _normalize_nucleotides(raw_nucleotides)

    if len(normalized_sequence) < 21:
        logger.error(f"Sequence too short: {len(normalized_sequence)} nt.")
        raise ValueError(
            f"Input sequence is extremely short ({len(normalized_sequence)} nt). "
            "A minimum length of 21 nucleotides is required to generate at least one viable siRNA candidate."
        )

    logger.info(f"Successfully loaded and normalized sequence of length {len(normalized_sequence)}.")
    return normalized_sequence
