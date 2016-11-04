"""
Provides viewer class for several pysofe objects.
"""

# IMPORTS
import numpy as np
from scipy.spatial import Delaunay

# backends
try:
    import matplotlib.pyplot as plt
except:
    pass

try:
    from . import canvas
except:
    pass

from .. import utils
from ..meshes.mesh import Mesh
from ..spaces.functions import FEFunction, MeshFunction

# DEBUGGING
from IPython import embed as IPS

BACKEND = 'vispy'

def use(backend):
    """
    Set visualisation backend.
    
    Parameters
    ----------

    backend : str ('vispy' | 'matplotlib')
        The visualisation backend to use.
    """

    # to set global module variable
    global BACKEND

    if backend.lower() in ('matplotlib', 'mpl', 'pyplot'):
        BACKEND = 'matplotlib'
    elif backend.lower() in ('vispy',):
        BACKEND = 'vispy'
    else:
        raise ValueError("Invalid visualisation backend ({})".format(backend))

def refgrid(dimension=2, resolution=1):
    """
    Creates a grid on the reference simplex of the given dimension.
    """

    # create equidistant points at lagrange nodes
    assert resolution > 0
    nodes = utils.lagrange_nodes(dimension=dimension,
                                 order=resolution)

    # triangulate these nodes
    cells = Delaunay(nodes.T).simplices.astype(np.uint32) + 1

    grid = Mesh(nodes=nodes.T, connectivity=cells)

    return grid

def visgrid(mesh, resolution=1):
    """
    Creates a visualisation grid as a refined version of
    the given mesh.

    Parameters
    ----------

    mesh : pysofe.meshes.mesh.Mesh
        The initial mesh

    resolution : int
        Local resolution of the visualization grid
    """

    # create reference grid
    rgrid = refgrid(dimension=mesh.dimension, resolution=resolution)

    # map reference grid points on given mesh
    vnodes = mesh.ref_map.eval(points=rgrid.nodes.T, deriv=0) # nE x nP x nD
    nE, nP, nD = vnodes.shape
    
    # create corresponding cell connectivity
    vcells = rgrid.cells[:,:,None] + (nP * np.arange(nE))
    
    vgrid = Mesh(nodes=np.vstack(vnodes),
                 connectivity=np.vstack(vcells.transpose([2,0,1])))

    return vgrid
    
class BaseViewerVP(object):
    """
    Base class for all viewers.
    """

    def __init__(self):
        self.fig = canvas.Figure(show=False)
        
    def _plot(self, *args, **kwargs):
        raise NotImplementedError()

    def show(self, *args, **kwargs):
        self._plot(*args, **kwargs)
        self.fig.show()

class MeshViewer(BaseViewerVP):
    """
    Viewer for :py:class:`pysofe.meshes.mesh.Mesh` instances.
    """

    def _plot(self, mesh, *args, **kwargs):
        # create axes
        ax = self.fig[0,0]

        if mesh.dimension is 1:
            x = mesh.nodes[:,0]
            y = np.zeros_like(x)

            ax.plot(x, y)
                
        elif mesh.dimension is 2:
            x, y = mesh.nodes.T

            ax.triplot(x, y, mesh.faces-1)
        elif mesh.dimension is 3:
            x, y, z = mesh.nodes.T

            ax.wireframe(x, y, z, mesh.faces-1)

        else:
            raise NotImplementedError('Invalid mesh dimension ({})'.format(mesh.dimension))

class DOFVectorViewer(BaseViewerVP):
    """
    Viewer for dof vectors.
    """

    def _plot(self, u, **kwargs):
        # get visualization data
        nodes, values, cells, edges = self._get_data(u, **kwargs)
        
        # setup axes
        n_values = values.shape[0]

        mode = kwargs.pop('mode', 'trisurface')
        edge_color = kwargs.get('edge_color', 'black')
        
        for i in xrange(n_values):
            ax = self.fig[0,i]

            if mode == 'trisurface':
                ax.trisurface(x=nodes[:,0], y=nodes[:,1], z=values[i],
                              triangles=cells - 1,
                              edges=edges - 1, edge_color=edge_color)
            elif mode == 'wireframe':
                ax.wireframe(x=nodes[:,0], y=nodes[:,1], z=values[i],
                             triangles=cells - 1,
                             edges=edges - 1)
            elif mode in ('triplot', 'heatmap'):
                raise NotImplementedError()

    def _get_data(self, u, **kwargs):
        # get visualization mesh
        local_res = kwargs.get('resolution', u.order+1)
        vmesh = visgrid(mesh=kwargs.get('mesh', u.fe_space.mesh),
                        resolution=local_res)
        
        # evaluate the dof vector
        deriv = kwargs.get('deriv', 0)
        
        if kwargs.get('eval_local', True):
            # create local reference grid
            rgrid = refgrid(dimension=vmesh.dimension, resolution=local_res)

            values = u(points=rgrid.nodes.T, deriv=deriv, local=True)

            if deriv is 0:
                values = values.ravel()
            elif deriv is 1:
                values = np.vstack(values).T
            else:
                raise ValueError("Invalid derivation order ({})".format(deriv))
        else:
            values = u(points=vmesh.nodes.T, deriv=deriv, local=False)
            values = values.T

        nodes = vmesh.nodes
        values = np.atleast_2d(values)
        cells = vmesh.cells
        edges = vmesh.edges.compress(vmesh.topology.get_boundary(d=1), axis=0)
        
        return nodes, values, cells, edges
        
class FunctionViewer(BaseViewerVP):
    """
    Viewer for callables.
    """

    def _plot(self, fnc, **kwargs):
        # get visualization data
        nodes, values, cells, edges = self._get_data(fnc, **kwargs)
        
        # setup axes
        n_values = values.shape[0]

        mode = kwargs.pop('mode', 'trisurface')
        
        for i in xrange(n_values):
            ax = self.fig[0,i]

            if mode == 'trisurface':
                ax.trisurface(x=nodes[:,0], y=nodes[:,1], z=values[i],
                              triangles=cells - 1,
                              edges=edges - 1, edge_color='black')
            elif mode == 'wireframe':
                ax.wireframe(x=nodes[:,0], y=nodes[:,1], z=values[i],
                             triangles=cells - 1,
                             edges=edges - 1)
            elif mode in ('triplot', 'heatmap'):
                pass

    def _get_data(self, fnc, **kwargs):
        # get visualization mesh
        local_res = kwargs.get('resolution', 1)

        if isinstance(fnc, MeshFunction):
            mesh = kwargs.get('mesh', fnc.mesh)
        else:
            mesh = kwargs.get('mesh', None)

            if mesh is None:
                raise RuntimeError("No mesh given!")

            fnc = MeshFunction(fnc, mesh)
        
        vmesh = visgrid(mesh=mesh, resolution=local_res)
        
        # evaluate the function on the visualization grid
        deriv = kwargs.get('deriv', 0)

        if deriv == 0:
            fnc_args = {}
        elif deriv == 1:
            # if not 'deriv' in self.fnc.fnc.func_code.co_varnames:
            #     raise RuntimeError("Function doesn't seem to support evaluation of derivatives!")
            # else:
            #     fnc_args = {'deriv', 1}
            fnc_args = {'deriv' : 1}
        else:
            raise ValueError("Invalid derivation order!")

        if kwargs.get('eval_local', True):
            # create local reference grid
            rgrid = refgrid(dimension=vmesh.dimension, resolution=local_res)

            # get flag specifying whether function is supposed to be continuous
            is_cont = kwargs.get('is_cont', True)
                
            if is_cont:
                values = fnc(points=rgrid.nodes.T, local=True, **fnc_args)
            else:
                # shift points towards barycenter of containing cell to ensure global
                # point search determines correct cell
                # (otherwise gps ambiguous for points on edges of physical mesh)
                bary = mesh.nodes.take(mesh.cells-1, axis=0).mean(axis=1).T
                
                # step size between points and barycenter of containing cell
                f = 1e-8
                
                ppoints = vmesh.nodes.T.reshape((vmesh.dimension, -1, rgrid.nodes.shape[0]))
                ppoints = (1-f) * ppoints + f * bary[:,:,None]
                ppoints = ppoints.reshape((vmesh.dimension, -1))
                
                # evaluate function in shifted points
                values = fnc(points=ppoints, local=False, **fnc_args)
                        
            if deriv == 0:
                values = values.ravel()
            elif deriv == 1:
                values = np.vstack(values).T
            else:
                raise ValueError("Invalid derivation order!")
        else:
            values = fnc(points=vmesh.nodes.T, local=False, **fnc_args)

        nodes = vmesh.nodes
        values = np.atleast_2d(values)
        cells = vmesh.cells
        edges = vmesh.edges.compress(vmesh.topology.get_boundary(d=1), axis=0)
        
        return nodes, values, cells, edges

