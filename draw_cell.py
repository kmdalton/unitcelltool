import numpy as np
import re
from numpy import sin,cos,sqrt
from pymol import cmd, cgo
from pymol.vfont import plain
from seaborn import color_palette

palette = color_palette('colorblind')

for i,rgb in enumerate(palette):
    cmd.set_color('color_{}'.format(chr(ord('a') + i)), rgb)

def get_orthogonalization_matrix(objectname):
    """
    get_orthogonalization_matrix(objectname : string)
    """
    a,b,c,alpha,beta,gamma,spacegroup = cmd.get_symmetry(objectname)
    alpha,beta,gamma=np.pi*alpha/180.,np.pi*beta/180.,np.pi*gamma/180.

    V = a*b*c*sqrt(1. - cos(alpha)*cos(alpha)
                      - cos(beta)*cos(beta)
                      - cos(gamma)*cos(gamma)
                      + 2*cos(alpha)*cos(beta)*cos(gamma))
    O = np.matrix([
        [a , b*cos(gamma), c*cos(beta)],
        [0., b*sin(gamma), c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)],
        [0.,           0., V/(a*b*sin(gamma))]
        ])

    return O

def draw_cell(objectname, length=10., colorname=None, aspect=None, origin=None):
    if aspect is not None:
        aspect = float(aspect)
    else:
        aspect = 0.03
    length = float(length) #For safety
    O = get_orthogonalization_matrix(objectname)
    T = np.matrix(cmd.get_object_matrix(objectname)).reshape((4,4))
    T1,R,T2 = T[3,:-1],T[:3,:3],T[:-1,3].T

    o = np.matrix([0., 0., 0.]).T
    a = np.matrix([1., 0., 0.]).T
    b = np.matrix([0., 1., 0.]).T
    c = np.matrix([0., 0., 1.]).T

    o = O*o
    a = O*a
    b = O*b
    c = O*c

    a = length*a/np.linalg.norm(a)
    b = length*b/np.linalg.norm(b)
    c = length*c/np.linalg.norm(c)

    o = T1 + o.T
    a = T1 + a.T
    b = T1 + b.T
    c = T1 + c.T

    o = R*o.T
    a = R*a.T
    b = R*b.T
    c = R*c.T

    o = T2 + o.T
    a = T2 + a.T
    b = T2 + b.T
    c = T2 + c.T


    if origin is not None:
        origin = re.sub(r'[^0-9,.]', '', origin)
        origin = np.array(map(float, origin.split(',')))
        a = a + (origin-o)
        b = b + (origin-o)
        c = c + (origin-o)
        o = origin

    a = list(np.array(a).flatten())
    b = list(np.array(b).flatten())
    c = list(np.array(c).flatten())
    o = list(np.array(o).flatten())

    radius = length*aspect

    obj = []
    if colorname is None:
        colorname = 'color_a'
    obj += cgo_arrow(o, a, aspect=aspect, color=colorname)
    obj += cgo_arrow(o, b, aspect=aspect, color=colorname)
    obj += cgo_arrow(o, c, aspect=aspect, color=colorname)
    cmd.load_cgo(obj, objectname + '_axis')
    cmd.pseudoatom(objectname + '_labels', label = 'a', color=cmd.get_color_tuple(colorname), pos = affine(a, o, 1.1))
    cmd.pseudoatom(objectname + '_labels', label = 'b', color=cmd.get_color_tuple(colorname), pos = affine(b, o, 1.1))
    cmd.pseudoatom(objectname + '_labels', label = 'c', color=cmd.get_color_tuple(colorname), pos = affine(c, o, 1.1))
    cmd.set("label_color", colorname, objectname + '_labels', )
    cmd.set("label_size", 30, objectname + '_labels')


def cgo_arrow(a, b, aspect, hlen=None, color=None):
    if hlen is None:
        hlen = 0.1
    if color is None:
        color = 'color_a'

    m = affine(a, b, hlen)
    length = np.linalg.norm(np.array(a) - np.array(b))
    color = cmd.get_color_tuple(color)
    color = list(color)
    obj = [cgo.CYLINDER] + a + m + [aspect*length] + color + color +\
          [cgo.CONE] + m + b + [0.5*hlen*length, 0.0] + color + color +\
          [1.0, 0.]

    return obj

def affine(a, b, theta):
    c = theta*np.array(a) + (1. - theta)*np.array(b)
    return list(np.array(c).flatten())

cmd.extend('draw_cell', draw_cell)
