# Advanced Usage: SPARC-X-API as a Socket Interface


## Overview
Experienced users can harness the power of SPARC and SPARC-X-API's
socket communication layer to build efficient and flexible
computational workflows. By integrating a socket communication
interface directly into SPARC, users can significantly reduce the
overhead typically associated with file I/O during calculation
restarts. This feature is particularly beneficial for tasks involving
repetitive operations like structural optimization and saddle point
searches, where traditional file-based communication can become a
bottleneck. The underlying software architecture is shown in
[Fig. 1](#scheme-sparc-socket):

```{figure} img/scheme_socket_hetero.png
:alt: scheme-sparc-socket
:name: scheme-sparc-socket

Fig. 1. SPARC electronic calculations with socket communication across hybrid computing platforms.
```

**Requirements**: the SPARC binary must be manually compiled from the source
code with socket
support and with the
`USE_SOCKET=1` flag enabled (see the [installation
instructions](https://github.com/SPARC-X/SPARC?tab=readme-ov-file#2-installation).

```{note}
You need SPARC C/C++ after version 2024.11.18 to enable socket support.
```

## Usage
The socket communication layer in SPARC and SPARC-X-API are designed for:
- **Efficiency:** Eliminates the need for intermediate file I/O, directly streaming data between processes.
- **Speed:** Enhances the performance of iterative calculations, crucial for large-scale simulations.
- **Flexibility:** Allows dynamic modification of calculation parameters without the need to restart the process.

The communication protocol implemented in SPARC and SPARC-X-API
adheres to the [i-PI protocol](https://github.com/i-pi/i-pi)
standard. Specifically, we implement the original i-PI protocol within
the SPARC C-source code, while the python SPARC-X-API uses a
backward-compatible protocol based on i-PI. The dual-mode design is
aimed for both low-level and high-level interfacing of the DFT codes,
providing the following features as shown in [Fig. 2](#scheme-sparc-protocol):

```{figure} img/scheme_sparc_protocol.png
:alt: scheme-sparc-protocol
:name: scheme-sparc-protocol

Fig. 2. Overview of the SPARC protocol as an extension to the standard i-PI protocol.
```

Based on the scenarios, the socket communication layer can be accessed
via the following approaches as shown in
[Fig. 3](#scheme-sparc-modes):

```{figure} img/scheme-SPARC-socket-modes.png
:name: scheme-sparc-modes
:alt: scheme-sparc-modes

Fig. 3. Different ways of using SPARC's socket mode.
```


1. **SPARC binary only** ([Fig. 3](#scheme-sparc-modes) **a**)

   SPARC binary with socket support can be readily coupled with any i-PI compatible socker server, such as
   `ase.calculators.socketio.SocketIOCalculator`, for example

   ```python
   from ase.calculators.socketio import SocketIOCalculator
   from subprocess import Popen
   calc = SocketIOCalculator(port=31415)
   with calc:
       # Start sparc at background
	   process = Popen("mpirun -n 8 sparc -name SPARC -socket localhost:31415", shell=True)
	   # Single point calculations
	   process.kill()
   ```

   The end user is responsible for generating the input files and
   making sure the same atoms structures are used by
   `SocketIOCalculator` and the SPARC binary. The mode is also limited
   to be run on a single computer system.


2. **Local-only mode** ([Fig. 3](#scheme-sparc-modes) **b**)

   Ideal for standalone calculations, this mode simulates a conventional calculator while benefiting from socket-based efficiency.

   ```python
   with SPARC(use_socket=True, **normal_parameters) as calc:
       # Execute single-point calculations
   ```
   For most users we recommend using this mode when performing a calculation on a single HPC node.

3. **Client (Relay) mode** ([Fig. 3](#scheme-sparc-modes) **c**)

   In this mode, the `sparc.SPARC` calculator servers as a passive
   client which listens to a remote i-PI-compatible server. When
   messages are received on the client side, it relays the relevant
   message to a local SPARC binary and send results back through the
   socket pipe. The server side can either be a normal i-PI compatible
   server (such as `SocketIOCalculator`) or server-mode `sparc.SPARC` (see 4).

   Start the client by:
   ```python
   client = SPARC(use_socket=True,
	              socket_params=dict(host="host.address.com", port=31415))
   with client:
       client.run()
   ```

   Or via Command-Line:
   ```bash
   python -m sparc.client -s host:port
   ```

   Note: when running SPARC-X-API as a socket client, the atoms object
   can be ommitted (is the server also runs the SPARC protocol). When
   new atoms positions and parameters arrive, the client will
   automatically determine if it is necessary to restart the SPARC
   subprocess.

4. **Server mode** ([Fig. 3](#scheme-sparc-modes) **d**)

   Paired with the client mode in (3), SPARC-X-API can be run as a
   socket server, isolated from the node that performs the
   computation. This can be useful for highly-distributed
   computational workflows.

   On the server node, run:
   ```python
   server_calc = SPARC(use_socket=True, socket_params=dict(port=31415, server_only=True), **normal_parameters)
   with server_calc:
       # Execute single point calculations for atoms_1
	   # Execute single point calculations for atoms_2
   ```

   In this case, the server will opens `0.0.0.0:31415` for
   connection. Make sure your server is directly accessible from the
   clients and the port is not occupied. The socker server is capable
   of receiving `raw_results` directly from the clients, making it
   possible to access `server_calc.raw_results` without access to the
   file systems on the client side.


## (In-progress) Controlling SPARC routines from socket interface

As shown in [Fig. 2](#scheme-sparc-protocol),
the SPARC socket protocol designs allows bidirectional control of
internal SPARC routines. Local- or server-mode `sparc.SPARC`
calculators can communicate with the SPARC binary via functions like
`set_params`. This can be useful for applications like on-the-fly
force field learning, electron density fitting, setting up boundary
conditions etc. Applications will be updated in both SPARC and
SPARC-X-API repositories.
