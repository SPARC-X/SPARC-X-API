import sys

from packaging import version

from sparc.api import SparcAPI  # Replace with your actual module import


def compare_versions():
    default_version = SparcAPI().sparc_version
    new_version = SparcAPI(json_api="parameters.json").sparc_version
    print(f"Default api version: {default_version}, upstream version {new_version}")

    if version.parse(default_version) < version.parse(new_version):
        print(f"Version changed from {default_version} to {new_version}")
        sys.exit(1)
    else:
        print("No version change")
        sys.exit(0)


if __name__ == "__main__":
    compare_versions()
