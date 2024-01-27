"""A i-PI compatible socket protocol implemented in SPARC
"""
import os
import random
import socket
import string
import io
import pickle

import numpy as np
from ase.calculators.socketio import (
    IPIProtocol,
    SocketClient,
    SocketServer,
    actualunixsocketname,
    SocketClosed,
)

import hashlib


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
        """Send an object dumped into pickle 
        """
        # We can use the highese protocol since the
        # python requirement >= 3.8
        pkl_bytes = pickle.dumps(obj, protocol=5)
        nbytes = len(pkl_bytes)
        md5_checksum = hashlib.md5(pkl_bytes)
        checksum_digest, checksum_count = (md5_checksum.digest(),
                                           md5_checksum.digest_size)
        self.log("  pickle bytes to send: ", str(nbytes))
        self.send(nbytes, np.int32)
        self.log("  sending pickle object....")
        self.socket.sendall(pkl_bytes)
        self.log(" sending md5 sum of size: ", str(checksum_count))
        self.send(checksum_count, np.int32)
        self.log(" sending md5 sum..... ", str(checksum_count))
        self.socket.sendall(checksum_digest)
        return

    def recv_object(self):
        """Return a decoded file
        """
        nbytes = int(self.recv(1, np.int32))
        self.log(" Will receive pickle object with n-bytes: ", nbytes)
        bytes_received = self._recvall(nbytes)
        checksum_nbytes = int(self.recv(1, np.int32))
        self.log(" Will receive cheksum digest of nbytes:", checksum_nbytes)
        digest_received = self._recvall(checksum_nbytes)
        digest_calc = hashlib.md5(bytes_received).digest()
        minlen = min(len(digest_calc), len(digest_received))
        assert digest_calc[:minlen] == digest_received[:minlen], ("MD5 checksum for the received object does not match!")
        obj = pickle.loads(bytes_received)
        return obj
    

    # def send_json(self, json_string, encoding="ascii"):
    #     """Send the full json-string in file mode.
    #     #TODO: add checksum
    #     """
    #     fd = io.BytesIO()
    #     # ASCII should be ok for current file system
    #     json_bytes = json_string.encode(encoding)
    #     with fd:
    #         self.log("Sending json-string in bytes mode")
    #         fd.seek(0)
    #         fd.write(json_bytes)
    #         fd.seek(0)
            
            
        

    def send_param(self, name, value):
        """Send a specific param setting to SPARC
        This is just a test function to see how things may work

        TODO:
        1) test with just 2 string values to see if SPARC can receive
        """
        self.log(f"Setup param {name}, {value}")
        msg = self.status()
        # TODO: see how NEEDINIT works
        # if msg == 'NEEDINIT':
        # self.sendinit()
        # msg = self.status()
        assert msg == "READY", msg
        # Send message
        self.sendmsg("SETPARAM")
        # Send name
        self.send_string(str(name))
        # Send value
        self.send_string(str(value))
        # After this step, socket client should return READY
        return

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
        self.parent_calc = parent_calc # Track the actual calculator
        # TODO: make sure the client is compatible with the default socketclient
        
        # We shall make NEEDINIT to be the default state
        self.state = "NEEDINIT"
        

    def calculate(self, atoms, use_stress):
        """Use the calculator instance
        """
        if atoms.calc is None:
            atoms.calc = self.parent_calc
        return super().calculate(atoms, use_stress)

    def irun(self, atoms, use_stress=True):
        """Reimplement single step calculation

        We're free to implement the INIT method in socket protocol as most
        calculators do not involve using these. We can let the C-SPARC to spit out
        error about needinit error.
        """
        try:
            while True:
                try:
                    msg = self.protocol.recvmsg()
                except SocketClosed:
                    # Server closed the connection, but we want to
                    # exit gracefully anyway
                    msg = 'EXIT'

                if msg == 'EXIT':
                    # Send stop signal to clients:
                    self.comm.broadcast(np.ones(1, bool), 0)
                    # (When otherwise exiting, things crashed and we should
                    # let MPI_ABORT take care of the mess instead of trying
                    # to synchronize the exit)
                    return
                elif msg == 'STATUS':
                    self.protocol.sendmsg(self.state)
                elif msg == 'POSDATA':
                    assert self.state == 'READY'
                    cell, icell, positions = self.protocol.recvposdata()
                    atoms.cell[:] = cell
                    atoms.positions[:] = positions

                    # User may wish to do something with the atoms object now.
                    # Should we provide option to yield here?
                    #
                    # (In that case we should MPI-synchronize *before*
                    #  whereas now we do it after.)

                    # Send signal for other ranks to proceed with calculation:
                    self.comm.broadcast(np.zeros(1, bool), 0)
                    energy, forces, virial = self.calculate(atoms, use_stress)

                    self.state = 'HAVEDATA'
                    yield
                elif msg == 'GETFORCE':
                    assert self.state == 'HAVEDATA', self.state
                    self.protocol.sendforce(energy, forces, virial)
                    self.state = 'NEEDINIT'
                elif msg == 'INIT':
                    assert self.state == 'NEEDINIT'
                    # At this step, we can ask the SPARC socket to 
                    # bead_index, initbytes = self.protocol.recvinit()
                    # self.bead_index = bead_index
                    # self.bead_initbytes = initbytes
                    self.state = 'READY'
                else:
                    raise KeyError('Bad message', msg)
        finally:
            self.close()
            
        def run(self, atoms=None, use_stress=False):
            """Socket mode in SPARC should allow arbitrary start
            """
            for _ in self.irun(atoms=atoms, use_stress=use_stress):
                pass
