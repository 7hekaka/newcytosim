#!/usr/bin/env python3
#
# FJN, Sainsbury Laboratory, Cambridge University
# 1.4.2021

"""
    Invoke 'play image' varying the rotation
    
Syntax:
    
    multi_views.py
    
Assemble images into a movie:

ffmpeg -i image%04d.tga -pix_fmt yuv420p -y multiviews.mp4

"""

import sys, math, os, subprocess


def quaternion(phi):
    """ Return Quaternion representing a rotation of angle phi"""
    # (theta, psi) define the axis of rotation:
    # to rotate around X, set q2=0 and q3=0
    # to rotate around Y, set q1=0 and q3=0
    # to rotate around Z, set q1=0 and q2=0
    theta = 0.5 * math.pi
    psi = 0.5 * math.pi
    q0 = math.cos(phi)
    q1 = math.sin(phi)*math.cos(theta)
    q2 = math.sin(phi)*math.sin(theta)*math.cos(psi)
    q3 = math.sin(phi)*math.sin(theta)*math.sin(psi)
    return [q0, q1, q2, q3]


def quaternion_repr(Q):
    """representation of a Quaternion"""
    s = ''
    c = '"'
    for x in Q:
        s += c
        s += repr(x)
        c = ' '
    return s+'"'


def make(inx, quat):
    """Generate one image"""
    exe = 'bin/play'
    rot = 'rotation='+quaternion_repr(quat)
    img = 'image%04i.tga' % inx
    sub = subprocess.Popen([exe, img, rot, 'axes=3'], stdout=subprocess.PIPE, text=True)
    out = ''
    while True:
        str = sub.stdout.readline()
        if not str:
            break;
        out = str
    sub.stdout.close()
    print(rot, out, end='')

#-------------------------------------------------------------------------------

def main(args):
    # process all orientations
    cnt = 0
    for i in range(0, 180, 10):
        make(cnt, quaternion(i*math.pi/180))
        cnt += 1


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])

