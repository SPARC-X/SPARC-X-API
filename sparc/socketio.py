"""A i-PI compatible socket protocol implemented in SPARC
"""
import numpy as np
import os
import socket
from ase.calculators.socketio import IPIProtocol, SocketServer, SocketClient
from ase.calculators.socketio import actualunixsocketname


class SPARCProtocol(IPIProtocol):
    """Extending the i-PI protocol to support extra routines
    """
    def send_string(self, msg, msglen=None):
        self.log('  send string', repr(msg))
        # assert msg in self.statements, msg
        if msglen is None:
            msglen = len(msg)
        assert msglen >= len(msg)
        msg = msg.encode('ascii').ljust(msglen)
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
        assert msg == 'READY', msg
        # Send message
        self.sendmsg("SETPARAM")
        # Send name
        self.send_string(str(name))
        # Send value
        self.send_string(str(value))
        # After this step, socket client should return READY
        return

    
class SPARCSocketServer(SocketServer):
    """We only implement the unix socket version due to simplicity"""
    #TODO: guard cases for non-unix sockets
    @property
    def socket_filename(self):
        return self.serversocket.getsockname()


class SPARCSocketClient(SocketClient):
    pass
