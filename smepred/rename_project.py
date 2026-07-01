import os
import re

ROOT_DIR = r"D:\Helixx\smepred"
EXCLUDE_DIRS = {'.git', 'models', '__pycache__', 'env', 'venv', '.venv'}
EXCLUDE_EXTS = {'.pkl', '.joblib', '.h5', '.pdf', '.png', '.jpg', '.zip', '.tar', '.gz'}

def replace_in_file(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError:
        return False

    original = content

    # Temporarily hide the GitHub URL to protect it
    placeholder_url = "GITHUB_REPO_URL_PLACEHOLDER"
    url_target = "github.com/nitinjadhav888/SMEpred-CDAC"
    content = content.replace(url_target, placeholder_url)

    # Perform replacements
    content = content.replace("SMEpred", "HelixZero-CMS")
    content = content.replace("SMEPred", "HelixZero-CMS")
    # For lowercase variables/imports like "smepred_existing", change to helixzero_cms
    content = content.replace("smepred", "helixzero_cms")

    # Restore the GitHub URL
    content = content.replace(placeholder_url, url_target)

    if content != original:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    return False

modified_files = []
for root, dirs, files in os.walk(ROOT_DIR):
    # Prune excluded directories
    dirs[:] = [d for d in dirs if d not in EXCLUDE_DIRS]
    
    for file in files:
        ext = os.path.splitext(file)[1].lower()
        if ext in EXCLUDE_EXTS:
            continue
        
        filepath = os.path.join(root, file)
        # Avoid modifying this script itself
        if file == "rename_project.py":
            continue
            
        if replace_in_file(filepath):
            modified_files.append(filepath)

print(f"Renamed occurrences in {len(modified_files)} files:")
for f in modified_files:
    print(f" - {os.path.relpath(f, ROOT_DIR)}")
