# IMPORTS
import numpy as np

from vispy import scene
from vispy.visuals import transforms

from . import visuals

class ViewWidget(scene.Widget):
    """
    Widget to facilitate drawing.

    Parameters
    ----------
    args : iterable
        Arguments passed to the `ViewBox` super class.

    kwargs : dict
        Keyword arguments passed to the `ViewBox` super class.
    """

    def __init__(self, *args, **kwargs):
        self.grid = None
        self.camera = None
        self.visuals = []
        self._cmap = get_colormap(name=kwargs.pop('cmap', 'hsl'))
        self._configured = False

        scene.Widget.__init__(self, *args, **kwargs)

        self.grid = self.add_grid(spacing=0, margin=10)
        
    def _configure(self, proj='3d'):
        if self._configured:
            return

        self.view = self.grid.add_view(row=0, col=0,)
                                       #bgcolor='#efefef',)
                                       #border_color='grey')

        if proj.lower() == '3d':
            self.view.camera = 'turntable'
            self.view.camera.set_range((-1,1), (-1,1), (-1,1))
        elif proj.lower() == '2d':
            self.view.camera = 'panzoom'
            self.view.camera.set_range((-1,1), (-1,1), (-1,1))
            self.view.camera.interactive = False
        else:
            raise ValueError("Invalid projection ({})".format(proj))
            
        self.camera = self.view.camera

        self._configured = True

    @property
    def cmap(self):
        """
        The colormap used for drawing.
        """
        return self._cmap

    @cmap.setter
    def cmap(self, colormap):
        self._cmap = get_colormap(colormap)
        self._update_vertex_colors()

    def _autoscale(self):
        if len(self.visuals) is 0:
            return
        
        abounds = np.array([self.axbounds(axis=i) for i in xrange(3)]).T

        a, b = abounds
        c, d = [-1., 1.]
        b_a = np.where(b - a == 0., 1., b - a)

        T = transforms.MatrixTransform()
        T.scale( (d - c) / b_a )
        T.translate( (b*c - a*d) / b_a )
        
        for mesh in self.visuals:
            mesh.transform = T
            
    def _update_vertex_colors(self):
        # get min/max values
        zmin, zmax = self.axbounds(axis=2)

        for mesh in self.visuals:
            if not isinstance(mesh, visuals.WireframeMeshVisual):
                # get values that determine the vertex color
                z = mesh.mesh_data.get_vertices()[:,2]
                
                # map them to [0,1]
                if not np.allclose(zmin, zmax):
                    t = np.expand_dims( (z - zmin) / (zmax - zmin) , axis=1)
                else:
                    t = np.expand_dims( (z - 0.5*zmin), axis=1)
                    
                # calculate new vertex colors
                vc = self.cmap.map(t)
                
                mesh.mesh_data.set_vertex_colors(vc)
                mesh.mesh_data_changed()

    def axbounds(self, axis):
        # get bounds of all visuals
        bounds = np.array([v.bounds(axis=axis) for v in self.visuals])
        if bounds.size > 0:
            lb = bounds[:,0].min()
            ub = bounds[:,1].max()
        else:
            lb, ub = ("", "")

        return lb, ub
        
    def plot(self, x, y, **kwargs):
        """
        Draw a data series using lines and markers.
        """

        self._configure(proj='2d')

        line = scene.LinePlot(data=(x, y), connect='strip', **kwargs)

        self.visuals += [line]
        self.view.add(line)
        self.view.camera.set_range()

        return line

    def triplot(self, x, y, triangles, **kwargs):
        """
        Draw an unstructured 2-dimensional grid.
        """

        self._configure(proj='2d')

        return self.wireframe(x=x, y=y, z=np.zeros_like(x),
                              triangles=triangles, **kwargs)
        
    def trisurface(self, x, y, z, triangles, **kwargs):
        """
        Draw a surface consisting of triangular patches.

        Parameters
        ----------

        x, y, z : array_like
            The vertex coordinate values as 1D arrays
        
        triangles : array_like
            The triangle connectivity array

        **kwargs : dict
            Keyword arguments to pass to :class:`SurfaceMeshVisual`
        """

        self._configure(proj='3d')
        
        vertices = np.column_stack([x, y, z]).astype(np.float32)
        faces = triangles.astype(np.uint32)
        edges = kwargs.pop('edges', None)

        color = kwargs.pop('color', None)
        vertex_colors = kwargs.pop('vertex_colors', None)
        face_colors = kwargs.pop('face_colors', None)
        edge_color = kwargs.pop('edge_color', None)

        zmin = kwargs.pop('zmin', z.min())
        zmax = kwargs.pop('zmax', z.max())

        if all(c is None for c in [color, vertex_colors, face_colors]):
            if not np.allclose(zmin, zmax):
                t = (z - zmin) / (zmax - zmin)
            else:
                t = (z - 0.5*zmin)

            vertex_colors = self.cmap.map( np.expand_dims(t, axis=1) )

        mesh = visuals.SurfaceMesh(vertices=vertices, faces=faces, edges=edges,
                                   vertex_colors=vertex_colors, face_colors=face_colors,
                                   edge_color=edge_color, color=color)

        self.visuals += [mesh]
        self._autoscale()
        self.view.add(mesh)
        
        return mesh

    def wireframe(self, x, y, z, triangles, **kwargs):
        """
        Draw a wireframe outline of a mesh.

        Parameters
        ----------

        x, y, z : array_like
            The vertex coordinate values as 1D arrays
        
        triangles : array_like
            The triangle connectivity array

        **kwargs : dict
            Keyword arguments to pass to :class:`WireframeMeshVisual`
        """

        self._configure(proj='3d')

        vertices = np.column_stack([x, y, z]).astype(np.float32)
        faces = triangles.astype(np.uint32)
        edges = kwargs.pop('edges', None)

        if kwargs.has_key('edge_color'):
            edge_color = kwargs.pop('edge_color', 'black')
        elif kwargs.has_key('color'):
            edge_color = kwargs.pop('color', 'black')
        else:
            edge_color = 'black'

        mesh = visuals.WireframeMesh(vertices=vertices, faces=faces, edges=edges,
                                     edge_color=edge_color)

        self.visuals += [mesh]
        self._autoscale()
        self.view.add(mesh)

        return mesh

    def _get_mesh_data(self, x, y, z, triangles, **kwargs):
        pass

# some wrappers for functional visualization
def trisurface(x, y, z, triangles, **kwargs):
    """
    See :py:func:`ViewWidget.trisurface`.
    """

    # instanciate canvas and get ViewWidget axes
    fig = Figure(show=True)
    ax = fig[0, 0]

    # use ViewWidget's method
    return ax.trisurface(x=x, y=y, z=z, triangles=triangles, **kwargs)

def heatmap(x, y, z, triangles, **kwargs):
    """
    See :py:func:`ViewWidget.heatmap`.
    """

    # instanciate canvas and get ViewWidget axes
    fig = Figure(show=True)
    ax = fig[0, 0]

    # use ViewWidget's method
    return ax.heatmap(x=x, y=y, z=z, triangles=triangles, **kwargs)

def wireframe(x, y, z, triangles, **kwargs):
    """
    See :py:func:`ViewWidget.wireframe`.
    """

    # instanciate canvas and get ViewWidget axes
    fig = Figure(show=True)
    ax = fig[0, 0]

    # use ViewWidget's method
    return ax.wireframe(x=x, y=y, z=z, triangles=triangles, **kwargs)
