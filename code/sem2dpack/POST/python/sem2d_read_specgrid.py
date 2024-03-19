import numpy as np

def sem2d_read_specgrid(datadir):
    if datadir and datadir[-1] != '/':
        datadir += '/'

    # Read mesh information
    nelem, npgeo, ngnod, npoin, ngll = np.loadtxt(datadir + 'grid_sem2d.hdr', dtype=int, skiprows=1)

    # Read finite element grid
    coorg = np.loadtxt(datadir + 'MeshNodesCoord_sem2d.tab')
    knods = np.loadtxt(datadir + 'ElmtNodes_sem2d.tab')

    # Read spectral element grid
    with open(datadir + 'coord_sem2d.dat', 'rb') as fid:
        coord = np.fromfile(fid, dtype=np.single).reshape((npoin, 2))

    with open(datadir + 'ibool_sem2d.dat', 'rb') as fid:
        ibool = np.fromfile(fid, dtype=np.int32).reshape((ngll, ngll, nelem))

    # Read GLL information
    with open(datadir + 'gll_sem2d.tab', 'r') as fid:
        data = np.loadtxt(fid)
    x = data[:, 0]  # GLL nodes in reference segment [-1:1]
    w = data[:, 1]  # GLL quadrature weights
    h = data[:, 2:].T  # Derivatives of Lagrange polynomials at the GLL nodes

    grid = {
        'nelem': nelem,
        'npgeo': npgeo,
        'ngnod': ngnod,
        'npoin': npoin,
        'ngll': ngll,
        'coorg': coorg,
        'knods': knods,
        'coord': coord,
        'ibool': ibool,
        'x': x,
        'w': w,
        'h': h
    }

    return grid

