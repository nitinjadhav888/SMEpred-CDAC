"""
sirna_generator.py — 21-mer siRNA Candidate Generation Engine

Responsible for parsing full mRNA transcripts and generating overlapping 
21-mer small interfering RNA (siRNA) candidate duplexes. 

Why 21-mers?
The RNA-induced silencing complex (RISC) specifically requires 21-nucleotide 
double-stranded RNAs to function efficiently. The engine physically slides a 
21-nt window across the target mRNA, generating the exact target (Sense) strand 
and deriving its reverse-complement (Antisense/Guide) strand for loading into Ago2.
"""

import logging
from dataclasses import dataclass
from typing import List, Dict, Any

logger = logging.getLogger(__name__)

# Translation table to compute the reverse complement of an RNA string
_RNA_COMPLEMENT = str.maketrans("AUGC", "UACG")


def _calculate_reverse_complement(sequence: str) -> str:
    """
    Computes the reverse-complement of an RNA sequence.
    
    Why: The Antisense (guide) strand is the exact reverse-complement of the Sense 
    (passenger) strand. The Ago2 protein loads the Antisense strand in a 5' to 3' 
    orientation to pair with the target mRNA. Therefore, we must complement the bases 
    and reverse the string to maintain the 5' -> 3' biological standard.
    
    Args:
        sequence (str): A 5' -> 3' RNA sequence (Sense strand).
        
    Returns:
        str: The 5' -> 3' reverse-complemented RNA sequence (Antisense strand).
    """
    return sequence.translate(_RNA_COMPLEMENT)[::-1]


@dataclass
class SiRNACandidate:
    """
    Represents a single 21-mer siRNA duplex candidate.
    
    Attributes:
        position (int): 0-based start index of the candidate within the parent mRNA.
        sense (str): The 5' -> 3' sequence perfectly matching the mRNA target region.
        antisense (str): The 5' -> 3' guide strand that will be loaded into RISC.
    """
    position: int
    sense: str
    antisense: str

    def to_dict(self) -> Dict[str, Any]:
        """Serializes the candidate to a JSON-compatible dictionary for API responses."""
        return {
            "position": self.position,
            "sense": self.sense,
            "antisense": self.antisense,
        }


def generate_candidates(mrna_sequence: str) -> List[SiRNACandidate]:
    """
    Generates an exhaustive list of all possible 21-mer siRNA candidates.
    
    Why: To find the optimal siRNA, we must evaluate every possible binding site on 
    the target mRNA. This function acts as a sliding window, moving 1 nucleotide at 
    a time to generate a complete combinatorial set of candidates for downstream 
    biophysical filtering and ML prediction.
    
    Args:
        mrna_sequence (str): The full, normalized mRNA/gene sequence.
        
    Returns:
        List[SiRNACandidate]: A complete list of all valid 21-mer siRNA pairs.
        
    Raises:
        ValueError: If the mRNA sequence is shorter than 21 nucleotides.
    """
    sirna_length = 21
    
    if len(mrna_sequence) < sirna_length:
        logger.error(f"Provided mRNA sequence is too short ({len(mrna_sequence)} nt).")
        raise ValueError(f"mRNA must be at least {sirna_length} nucleotides long.")

    candidates: List[SiRNACandidate] = []
    total_candidates = len(mrna_sequence) - sirna_length + 1
    
    for i in range(total_candidates):
        sense_strand = mrna_sequence[i : i + sirna_length]
        antisense_strand = _calculate_reverse_complement(sense_strand)
        candidates.append(
            SiRNACandidate(
                position=i, 
                sense=sense_strand, 
                antisense=antisense_strand
            )
        )

    logger.info(f"Generated {len(candidates)} raw 21-mer candidates from mRNA input.")
    return candidates


def generate_dsirna_candidate(dsirna_sequence: str) -> List[SiRNACandidate]:
    """
    Extracts the single active 21-mer from a 25–30 nt Dicer-substrate RNA (DsiRNA).
    
    Why: Dicer is an endogenous enzyme that processes longer double-stranded RNAs. 
    It anchors at the 3' end and cleaves exactly ~21 nucleotides to produce a mature 
    siRNA. This function mimics Dicer cleavage by extracting the terminal 21-mer 
    from a user-provided DsiRNA sequence, allowing the model to predict the efficacy 
    of the final biological product.
    
    Args:
        dsirna_sequence (str): The DsiRNA sequence (25–30 nt).
        
    Returns:
        List[SiRNACandidate]: A single-element list containing the mature 21-mer product.
        
    Raises:
        ValueError: If the input is not within the biological 25-30 nt DsiRNA range.
    """
    sirna_length = 21
    
    if not (25 <= len(dsirna_sequence) <= 30):
        logger.error(f"Invalid DsiRNA length: {len(dsirna_sequence)} nt.")
        raise ValueError(f"DsiRNA input must be 25–30 nt, got {len(dsirna_sequence)}.")
        
    # Mimic Dicer cleavage: Extract the 5' 21-mer
    sense_strand = dsirna_sequence[:sirna_length]
    antisense_strand = _calculate_reverse_complement(sense_strand)
    
    logger.info(f"Successfully processed DsiRNA sequence. Extracted mature 21-mer.")
    return [SiRNACandidate(position=0, sense=sense_strand, antisense=antisense_strand)]
