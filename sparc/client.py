import argparse
from pathlib import Path

from ase.io import read

import sparc

from .calculator import SPARC
from .io import read_sparc


def split_socket_name(name):
    """Split the host:address

    Returns:
      (host, port, unixsocket) tuple
    """
    name = name.strip()
    assert ":" in name, "Socket name must be in the format of host:port"
    host, port = name.split(":")
    if port.lower() == "unix":
        unixsocket = host
        port = -1
        host = ""
        assert len(unixsocket) > 0, "Unix socket length must be non-zero!"
    else:
        port = int(port)
        unixsocket = None
        host = host if len(host) > 0 else "localhost"

    return host, port, unixsocket


def main():
    """Running SPARC-X-API as a socket client from command line

    The client wraps a socket communication layer on top of the SPARC
    binary. The implementation is both compatible with
    `SocketIOCalculator` in ase (using standard i-PI protocol), or
    with SPARC-X-API in socket server mode (using SPARC's extended
    i-PI protocol).

    Usage:
    1. Start the socket client and outputs to current directory:
        python -m sparc.client -s host:port --workdir .

       If the workdir is a SPARC calculation bundle, the initial atoms and parameters will be reloaded.

    2. Start the socket client with initial atoms read from file
        python -m sparc.client -s host:port --atoms-from-file input.xyz

    If the client is communicating with the standard i-PI server, an
    initial atoms object is required (either read from SPARC input
    files or via --atoms-from-file). However, if the server uses the
    SPARC protocol, the client can be started without initial atoms.
    """
    parser = argparse.ArgumentParser(
        usage=main.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--socket",
        help=(
            "Address of the socket server "
            "in the format of host:port. If host is not defined, localhost will be used."
        ),
    )
    parser.add_argument(
        "-w",
        "--workdir",
        help=("Workdir for performing the SPARC calculations"),
        default=".",
    )
    parser.add_argument(
        "--atoms-from-file",
        help=("File or directory to read the input atoms information"),
        default=None,
    )
    parser.add_argument(
        "--atoms-format",
        help="File format to read from external file.",
        default=None,
    )
    args = parser.parse_args()
    host, port, unixsocket = split_socket_name(args.socket)
    workdir = Path(args.workdir)
    # TODO: implement unixsocket
    # TODO: reuse init params
    try:
        init_atoms = read_sparc(workdir)
    except:
        init_atoms = None
    if (init_atoms is None) and args.atoms_from_file:
        atoms_file = Path(args.atoms_from_file)
        atoms_format = args.atoms_format
        init_atoms = read(atoms_file, format=atoms_format)

    client_calc = SPARC(
        directory=workdir,
        use_socket=True,
        socket_params=dict(host=host, port=port, server_only=False),
    )
    # We should always enable use_stress, since the
    # SPARC implementation ensures it will also be present
    client_calc.run_client(atoms=init_atoms, use_stress=True)


if __name__ == "__main__":
    main()
