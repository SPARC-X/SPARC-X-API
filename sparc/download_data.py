"""Download the pseudopotential and other related files after sparc-x-api is installed

Run:

python -m sparc.download_data
"""

import hashlib
import shutil
import tempfile
import zipfile
from io import BytesIO
from pathlib import Path

# import urllib.request
from urllib.request import urlopen

from .common import psp_dir

sparc_tag = "b702c1061400a2d23c0e223e32182609d7958156"
sparc_source_url = "https://github.com/SPARC-X/SPARC/archive/{sparc_tag}.zip"
# This is a all_psp8_checksum
all_psp8_checksum = "5ef42c4a81733a90b0e080b771c5a73a"


def download_psp(sparc_tag=sparc_tag, psp_dir=psp_dir):
    """Download the external PSPs into the sparc/psp folder

    Arguments:
        sparc_tag (str): Commit hash or git tag for the psp files
        psp_dir (str or PosixPath): Directory to download the psp files
    """
    if is_psp_download_complete():
        print("PSPs have been successfully downloaded!")
        return
    download_url = sparc_source_url.format(sparc_tag=sparc_tag)
    print(f"Download link: {download_url}")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        with urlopen(download_url) as zipresp:
            with zipfile.ZipFile(BytesIO(zipresp.read())) as zfile:
                zfile.extractall(tmpdir)
                # print(list(os.walk(tmpdir)))
                source_dir = next(tmpdir.glob("SPARC-*/psps"))
                print(f"Found source_dir at {source_dir}")
                if not source_dir.is_dir():
                    raise FileNotFoundError("Error downloading or extracting zip")
                print(f"Moving psp files to {psp_dir}")
                for ext in ("*.psp8", "*.psp", "*.pot"):
                    for pspf in source_dir.glob(ext):
                        print(f"Found {pspf} --> {psp_dir}")
                        shutil.copy(pspf, psp_dir)
    if not is_psp_download_complete(psp_dir):
        raise RuntimeError(f"Files downloaded to {psp_dir} have different checksums!")
    return


def checksum_all(psp_dir=psp_dir, extension="*.psp8"):
    """Checksum all the files under the psp_dir to make sure the psp8 files
    are the same as intended

    Arguments:
        psp_dir (str or PosixPath): Directory for the psp files
        extension (str): Search pattern for the psp files, either '.psp', '.psp8' or '.pot'

    Returns:
        str: Checksum for all the files concatenated
    """
    checker = hashlib.md5()
    psp_dir = Path(psp_dir)
    # Use sorted to make sure file order is correct
    for filename in sorted(psp_dir.glob(extension)):
        # Open the file in binary mode and update the group checksum
        with open(filename, "r") as f:
            f_checker = hashlib.md5()
            content = f.read().encode("utf8")
            f_checker.update(content)
            checker.update(f_checker.hexdigest().encode("ascii"))
    final_checksum = checker.hexdigest()
    return final_checksum


def is_psp_download_complete(psp_dir=psp_dir):
    return checksum_all(psp_dir) == all_psp8_checksum


if __name__ == "__main__":
    print("Running command-line psp downloader")
    download_psp()
