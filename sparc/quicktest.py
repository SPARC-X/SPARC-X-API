"""A simple test module for SPARC
"""

def main():
    from .calculator import SPARC
    from .sparc_io_bundle import SparcBundle
    calc = SPARC()
    calc._make_command()
    print("ALL PASS!")
    
if __name__ == "__main__":
    main()