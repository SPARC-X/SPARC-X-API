"""A i-PI compatible socket protocol implemented in SPARC
"""
import hashlib
import io
import os
import pickle
import random
import socket
import string

import numpy as np
from ase.calculators.socketio import (
    IPIProtocol,
    SocketClient,
    SocketClosed,
    SocketServer,
    actualunixsocketname,
)


def generate_random_socket_name(prefix="sparc_", length=6):
    """Generate a random socket name with the given prefix and a specified length of random hex characters."""
    random_chars = "".join(random.choices(string.hexdigits.lower(), k=length))
    return prefix + random_chars


class SPARCProtocol(IPIProtocol):
    """Extending the i-PI protocol to support extra routines"""

    def send_string(self, msg, msglen=None):
        self.log("  send string", repr(msg))
        # assert msg in self.statements, msg
        if msglen is None:
            msglen = len(msg)
        assert msglen >= len(msg)
        msg = msg.encode("ascii").ljust(msglen)
        self.send(msglen, np.int32)
        self.socket.sendall(msg)
        return

    def send_object(self, obj):
        """Send an object dumped into pickle"""
        # We can use the highese protocol since the
        # python requirement >= 3.8
        pkl_bytes = pickle.dumps(obj, protocol=5)
        nbytes = len(pkl_bytes)
        md5_checksum = hashlib.md5(pkl_bytes)
        checksum_digest, checksum_count = (
            md5_checksum.digest(),
            md5_checksum.digest_size,
        )
        self.sendmsg("PKLOBJ")  # To distinguish from other methods like INIT
        self.log("  pickle bytes to send: ", str(nbytes))
        self.send(nbytes, np.int32)
        self.log("  sending pickle object....")
        self.socket.sendall(pkl_bytes)
        self.log(" sending md5 sum of size: ", str(checksum_count))
        self.send(checksum_count, np.int32)
        self.log(" sending md5 sum..... ", str(checksum_count))
        self.socket.sendall(checksum_digest)
        return

    def recv_object(self, include_header=True):
        """Return a decoded file

        include_header: should we receive the header or not
        """
        if include_header:
            msg = self.recvmsg()
            assert (
                msg.strip() == "PKLOBJ"
            ), f"Incorrect header {msg} received when calling recv_object method! Please contact the developers"
        nbytes = int(self.recv(1, np.int32))
        self.log(" Will receive pickle object with n-bytes: ", nbytes)
        bytes_received = self._recvall(nbytes)
        checksum_nbytes = int(self.recv(1, np.int32))
        self.log(" Will receive cheksum digest of nbytes:", checksum_nbytes)
        digest_received = self._recvall(checksum_nbytes)
        digest_calc = hashlib.md5(bytes_received).digest()
        minlen = min(len(digest_calc), len(digest_received))
        assert (
            digest_calc[:minlen] == digest_received[:minlen]
        ), "MD5 checksum for the received object does not match!"
        obj = pickle.loads(bytes_received)
        return obj

    def send_param(self, name, value):
        """Send a specific param setting to SPARC
        This is just a test function to see how things may work

        TODO:
        1) test with just 2 string values to see if SPARC can receive
        """
        self.log(f"Setup param {name}, {value}")
        msg = self.status()
        assert msg == "READY", msg
        # Send message
        self.sendmsg("SETPARAM")
        # Send name
        self.send_string(str(name))
        # Send value
        self.send_string(str(value))
        # After this step, socket client should return READY
        return

    def sendinit(self):
        """Mimick the old sendinit method but to provide atoms and params
        to the calculator instance.
        The actual behavior regarding how the calculator would be (re)-initialized, dependends on the implementation of recvinit
        """
        self.log(" New sendinit for SPARC protocol")
        self.sendmsg("INIT")
        self.send(0, np.int32)  # fallback
        msg_chars = [ord(c) for c in "NEWPROTO"]
        len_msg = len(msg_chars)
        self.send(len_msg, np.int32)
        self.send(msg_chars, np.byte)  # initialization string
        return

    def recvinit(self):
        """Fallback recvinit method"""
        return super().recvinit()

    def calculate_new_protocol(self, atoms, params):
        atoms = atoms.copy()
        atoms.calc = None
        self.log(" calculate with new protocol")
        msg = self.status()
        # We don't know how NEEDINIT is supposed to work, but some codes
        # seem to be okay if we skip it and send the positions instead.
        if msg == "NEEDINIT":
            self.sendinit()
            self.send_object((atoms, params))
            msg = self.status()
        cell = atoms.get_cell()
        positions = atoms.get_positions()  # Original order
        assert msg == "READY", msg
        icell = np.linalg.pinv(cell).transpose()
        self.sendposdata(cell, icell, positions)
        msg = self.status()
        assert msg == "HAVEDATA", msg
        e, forces, virial, morebytes = self.sendrecv_force()
        r = dict(energy=e, forces=forces, virial=virial, morebytes=morebytes)
        # Additional data (e.g. parsed from file output)
        moredata = self.recv_object()
        return r, moredata


# TODO: make sure both calc are ok


class SPARCSocketServer(SocketServer):
    """We only implement the unix socket version due to simplicity

    parent: the SPARC parent calculator
    """

    def __init__(
        self,
        port=None,
        unixsocket=None,
        timeout=None,
        log=None,
        parent=None
        # launch_client=None,
    ):
        super().__init__(port=port, unixsocket=unixsocket, timeout=timeout, log=log)
        self.parent = parent
        print("Parent : ", self.parent)
        if self.parent is not None:
            self.proc = self.parent.process
        else:
            self.proc = None
        print(self.proc)

    # TODO: guard cases for non-unix sockets
    @property
    def socket_filename(self):
        return self.serversocket.getsockname()

    @property
    def proc(self):
        if self.parent:
            return self.parent.process
        else:
            return None

    @proc.setter
    def proc(self, value):
        return

    def _accept(self):
        """Use the SPARCProtocol instead"""
        print(self.proc)
        super()._accept()
        print(self.proc)
        old_protocol = self.protocol
        # Swap the protocol
        if old_protocol:
            self.protocol = SPARCProtocol(self.clientsocket, txt=self.log)
        return

    def send_atoms_and_params(self, atoms, params):
        """Update the atoms and parameters for the SPARC calculator
        The params should be assignable to SPARC.set

        The calc for atoms is stripped for simplicity
        """
        atoms.calc = None
        params = dict(params)
        pair = (atoms, params)
        self.protocol.send_object(pair)
        return

    def calculate_origin_protocol(self, atoms):
        """Send geometry to client and return calculated things as dict.

        This will block until client has established connection, then
        wait for the client to finish the calculation."""
        assert not self._closed

        # If we have not established connection yet, we must block
        # until the client catches up:
        if self.protocol is None:
            self._accept()
        return self.protocol.calculate(atoms.positions, atoms.cell)

    def calculate_new_protocol(self, atoms, params={}):
        assert not self._closed

        # If we have not established connection yet, we must block
        # until the client catches up:
        if self.protocol is None:
            self._accept()
        return self.protocol.calculate_new_protocol(atoms, params)


class SPARCSocketClient(SocketClient):
    def __init__(
        self,
        host="localhost",
        port=None,
        unixsocket=None,
        timeout=None,
        log=None,
        parent_calc=None
        # use_v2_protocol=True    # If we should use the v2 SPARC protocol
    ):
        """Reload the socket client and use SPARCProtocol"""
        super().__init__(
            host=host,
            port=port,
            unixsocket=unixsocket,
            timeout=timeout,
            log=log,
        )
        sock = self.protocol.socket
        self.protocol = SPARCProtocol(sock, txt=log)
        self.parent_calc = parent_calc  # Track the actual calculator
        # TODO: make sure the client is compatible with the default socketclient

        # We shall make NEEDINIT to be the default state
        # self.state = "NEEDINIT"

    def calculate(self, atoms, use_stress):
        """Use the calculator instance"""
        if atoms.calc is None:
            atoms.calc = self.parent_calc
        return super().calculate(atoms, use_stress)

    def irun(self, atoms, use_stress=True):
        """Reimplement single step calculation

        We're free to implement the INIT method in socket protocol as most
        calculators do not involve using these. We can let the C-SPARC to spit out
        error about needinit error.
        """
        # Discard positions received from POSDATA
        # if the server has send positions through recvinit method
        discard_posdata = False
        new_protocol = False
        try:
            while True:
                try:
                    msg = self.protocol.recvmsg()
                except SocketClosed:
                    # Server closed the connection, but we want to
                    # exit gracefully anyway
                    msg = "EXIT"

                if msg == "EXIT":
                    # Send stop signal to clients:
                    self.comm.broadcast(np.ones(1, bool), 0)
                    # (When otherwise exiting, things crashed and we should
                    # let MPI_ABORT take care of the mess instead of trying
                    # to synchronize the exit)
                    return
                elif msg == "STATUS":
                    self.protocol.sendmsg(self.state)
                elif msg == "POSDATA":
                    assert self.state == "READY"
                    assert (
                        atoms is not None
                    ), "Your SPARCSocketClient isn't properly initialized!"
                    cell, icell, positions = self.protocol.recvposdata()
                    if not discard_posdata:
                        atoms.cell[:] = cell
                        atoms.positions[:] = positions

                    # At this stage, we should only rely on self.calculate
                    # to continue the socket calculation or restart
                    self.comm.broadcast(np.zeros(1, bool), 0)
                    energy, forces, virial = self.calculate(atoms, use_stress)

                    self.state = "HAVEDATA"
                    yield
                elif msg == "GETFORCE":
                    assert self.state == "HAVEDATA", self.state
                    self.protocol.sendforce(energy, forces, virial)
                    if new_protocol:
                        # TODO: implement more raw results
                        raw_results = self.parent_calc.raw_results
                        self.protocol.send_object(raw_results)
                    self.state = "NEEDINIT"
                elif msg == "INIT":
                    assert self.state == "NEEDINIT"
                    # Fall back to the default socketio
                    bead_index, initbytes = self.protocol.recvinit()
                    # The parts below use the new sparc protocol
                    print("Init bytes: ", initbytes)
                    init_msg = "".join([chr(d) for d in initbytes])
                    if init_msg.startswith("NEWPROTO"):
                        new_protocol = True
                        recv_atoms, params = self.protocol.recv_object()
                        print(recv_atoms, params)
                        if params != {}:
                            self.parent_calc.set(**params)
                        # TODO: should we update the atoms directly or keep copy?
                        atoms = recv_atoms
                        atoms.calc = self.parent_calc
                        discard_posdata = True
                    self.state = "READY"
                else:
                    raise KeyError("Bad message", msg)
        finally:
            self.close()

    def run(self, atoms=None, use_stress=False):
        """Socket mode in SPARC should allow arbitrary start"""
        # As a default we shall start the SPARCSocketIO always in needinit mode
        if atoms is None:
            self.state = "NEEDINIT"
        for _ in self.irun(atoms=atoms, use_stress=use_stress):
            pass
