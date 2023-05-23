def main(prog='sparc-ase',
         description='ASE command line tool with SPARC support',
         hook=None, args=None):
    import sparc
    from ase.cli.main import main
    main(prog=prog, description=description, hook=hook, args=args)
    

if __name__ == "__main__":
    main()
