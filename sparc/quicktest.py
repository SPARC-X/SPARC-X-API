"""A simple test module for SPARC
"""


def main():
    from sparc import SPARC

    calc = SPARC()
    calc._make_command()
    print("ALL PASS!")


if __name__ == "__main__":
    main()
