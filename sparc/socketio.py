"""A i-PI compatible socket protocol implemented in SPARC
"""
import os
import socket

import numpy as np
from ase.calculators.socketio import (
    IPIProtocol,
    SocketClient,
    SocketServer,
    actualunixsocketname,
)


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
        comm=None,
    ):
        """Reload the socket client and use SPARCProtocol"""
        super().__init__(
            host=host,
            port=port,
            unixsocket=unixsocket,
            timeout=timeout,
            log=log,
            comm=comm,
        )
        sock = self.protocol.socket
        self.protocol = SPARCProtocol(sock, txt=log)
