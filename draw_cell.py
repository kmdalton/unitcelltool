import numpy as np
from numpy import sin,cos,sqrt
from pymol import cmd, cgo, CmdException


colors = {
    'axis_a': [27/255.,158/255.,119/255.],
    'axis_b': [217/255.,95/255.,2/255.],
    'axis_c': [117/255.,112/255.,179/255.]
}

for k in colors:
    cmd.set_color(k, colors[k])

def get_orthogonalization_matrix(objectname):
    """
    get_orthogonalization_matrix(objectname -- string)
    """
    print cmd.get_symmetry(objectname)
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

def draw_cell(objectname, length=10.):
    length = float(length) #For safety
    O = get_orthogonalization_matrix(objectname)
    T = np.matrix(cmd.get_object_matrix(objectname)).reshape((4,4))
    T1,R,T2 = T[3,:-1],T[:3,:3],T[:-1,3].T
    print cmd.get_object_matrix(objectname)


    o = np.matrix([0., 0., 0.]).T
    a = np.matrix([1., 0., 0.]).T
    b = np.matrix([0., 1., 0.]).T
    c = np.matrix([0., 0., 1.]).T

    o = O*o
    a = O*a
    b = O*b
    c = O*c

    print a.shape
    a = length*a/np.linalg.norm(a)
    b = length*b/np.linalg.norm(b)
    c = length*c/np.linalg.norm(c)

    print a.shape

    o = T1 + o.T
    a = T1 + a.T
    b = T1 + b.T
    c = T1 + c.T

    print a.shape

    o = R*o.T
    a = R*a.T
    b = R*b.T
    c = R*c.T

    print a.shape

    o = T2 + o.T
    a = T2 + a.T
    b = T2 + b.T
    c = T2 + c.T

    print a.shape
    a = list(np.array(a).flatten())
    b = list(np.array(b).flatten())
    c = list(np.array(c).flatten())
    o = list(np.array(o).flatten())

    print a
    print b
    print c

    cgo_arrow(o, a, color='axis_a', name=objectname + "_a")
    cgo_arrow(o, b, color='axis_b', name=objectname + "_b")
    cgo_arrow(o, c, color='axis_c', name=objectname + "_c")

'''
http://pymolwiki.org/index.php/cgo_arrow

(c) 2013 Thomas Holder, Schrodinger Inc.

License: BSD-2-Clause
'''

def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
              color='blue red', name=''):
    '''
DESCRIPTION

    Create a CGO arrow between two picked atoms.

ARGUMENTS

    atom1 = string: single atom selection or list of 3 floats {default: pk1}

    atom2 = string: single atom selection or list of 3 floats {default: pk2}

    radius = float: arrow radius {default: 0.5}

    gap = float: gap between arrow tips and the two atoms {default: 0.0}

    hlength = float: length of head

    hradius = float: radius of head

    color = string: one or two color names {default: blue red}

    name = string: name of CGO object
    '''
    from chempy import cpv

    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)

    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
          [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
          [1.0, 0.0]

    if not name:
        name = cmd.get_unused_name('arrow')

    cmd.load_cgo(obj, name)

cmd.extend('cgo_arrow', cgo_arrow)
cmd.extend('draw_cell', draw_cell)
