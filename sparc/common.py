from pathlib import Path

import pkg_resources

repo_dir = Path(pkg_resources.resource_filename("sparc", "."))
psp_dir = Path(pkg_resources.resource_filename("sparc", "psp"))
