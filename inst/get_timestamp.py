#!/usr/bin/env python3

""" get timestamp
Get an NTP timestamp.
"""

import socket
import struct

def gettime_ntp(addr='time.nist.gov'):
    """ gettime_ntp
        Get NTP Timestamp for file versioning.
        Arguments
            * addr (str) : valid NTP address (e.g. '0.uk.pool.ntp.org',
                'time.nist.gov' etc)
        Returns
            * timestamp (str) : NTP seconds timestamp, converted to string.
    """
    TIME1970 = 2208988800
    client = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    data = '\x1b' + 47 * '\0'
    client.sendto(data.encode('utf-8'), (addr, 123))
    data, address = client.recvfrom(1024)
    t = struct.unpack('!12I', data)[10] - TIME1970
    return str(t)

if __name__ == "__main__":
    print("Getting timestamp...")
    gettime_ntp()