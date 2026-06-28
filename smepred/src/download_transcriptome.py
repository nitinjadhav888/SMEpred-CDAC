import urllib.request
import gzip
import shutil
import os
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

URL = "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "data"
GZ_FILE = DATA_DIR / "Homo_sapiens.GRCh38.cdna.all.fa.gz"
FASTA_FILE = DATA_DIR / "human_transcriptome.fasta"

def download_and_extract():
    DATA_DIR.mkdir(exist_ok=True)
    
    if FASTA_FILE.exists():
        logger.info(f"{FASTA_FILE} already exists. Skipping download.")
        return

    logger.info(f"Downloading {URL}...")
    try:
        urllib.request.urlretrieve(URL, GZ_FILE)
        logger.info("Download complete. Extracting...")
        
        with gzip.open(GZ_FILE, 'rb') as f_in:
            with open(FASTA_FILE, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                
        logger.info(f"Extraction complete: {FASTA_FILE}")
        
        # Clean up zip
        if GZ_FILE.exists():
            os.remove(GZ_FILE)
            
    except Exception as e:
        logger.error(f"Failed to download/extract: {e}")

if __name__ == "__main__":
    download_and_extract()
