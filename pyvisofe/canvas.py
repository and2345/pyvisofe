"""
Provides a canvas to draw on.
"""

# IMPORTS
import numpy as np

from vispy import scene
from vispy.visuals import transforms
from vispy.color import get_colormap

from . import visuals

# DEBUGGING
from IPython import embed as IPS

class Figure(scene.SceneCanvas):
    """
    Custom figure class.
    """

    def __init__(self, size=(800, 600), bgcolor='#efefef', show=False, **kwargs):
        # initialiaze before freeze occurs
        self._grid = None
        self._widgets = []

        self._cmap = get_colormap(name=kwargs.pop('cmap', 'hsl'))
        self._cbar = None

        scene.SceneCanvas.__init__(self, size=size, bgcolor=bgcolor, show=show,
                                   keys='interactive', **kwargs)

        self._grid = self.central_widget.add_grid()
        self._grid._default_class = ViewWidget

    def __getitem__(self, ids):
        vw = self._grid.__getitem__(ids)
        self._widgets += [vw]
        return vw

    @property
    def widgets(self):
        """
        List of associated ViewWidget instances.
        """
        return tuple(self._widgets)

    @property
    def cmap(self):
        """
        The colormap used for drawing.
        """
        return self._cmap

    @cmap.setter
    def cmap(self, colormap):
        self._cmap = get_colormap(colormap)

        if self._cbar is not None:
            self._cbar.cmap = self.cmap

        for widget in self._widgets:
            widget.cmap = self.cmap

    def colorbar(self):
        """
        Adds a colorbar to the figure.
        """
        
        if self._cbar is None:
            # create colorbar widget
            zbounds = np.array([widget.axbounds(axis=-1) for widget in self.widgets])
            clim = (zbounds[:,0].min(), zbounds[:,1].max())
            cbar = scene.ColorBarWidget(cmap=self.cmap, clim=clim,
                                        orientation='left')
            cbar.stretch = [0.2, 1.0]
            
            self._cbar = cbar
            self._grid.add_widget(cbar)
        else:
            pass
            
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
        
    def _configure(self, proj='3d', interactive=True):
        if self._configured:
            return

        self.view = self.grid.add_view(row=0, col=0,)
                                       #bgcolor='#efefef',)
                                       #border_color='grey')

        if proj.lower() == '3d':
            self.view.camera = 'turntable'
            self.view.camera.set_range((-1,1), (-1,1), (-1,1))
            self.view.camera.interactive = interactive
        elif proj.lower() == '2d':
            self.view.camera = 'panzoom'
            self.view.camera.set_range((-1,1), (-1,1), (-1,1))
            self.view.camera.interactive = interactive
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

    def _get_vertex_colors(self, values, vmin=None, vmax=None, alpha=1.):
        assert values.ndim == 1

        vmin = vmin if vmin is not None else values.min()
        vmax = vmax if vmax is not None else values.max()

        if not np.allclose(vmin, vmax):
            t = (values - vmin) / (vmax - vmin)
        else:
            t = (values - 0.5*vmin)
            
        vcolors = self.cmap.map( np.expand_dims(t, axis=1) )
        vcolors[:,3] = alpha

        return vcolors
            
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


    def plot(self, x, y, z=None, color='black', symbol='o',
             width=1., marker_size=10., connect='strip'):
        """
        Draw a data series using lines and markers.

        Parameters
        ----------

        x, y[, z] : array_like
            The data values

        color : str
            The color of the lines

        symbol : str
            Marker symbol to use

        width : float
            The line width in px

        marker_size : float
           The marker size in px (if `0` no markers will be shown)
    
        connect : ['strip'] | 'segments'
            Determines how the data points are connected by lines
        """
        self.lineplot(x, y, z, color, symbol, width, marker_size, connect)

    def lineplot(self, x, y, z=None, connect='strip',
                 color='black', symbol='o', width=1., marker_size=10.):
        """
        Draw a series of lines and markers.

        Parameters
        ----------

        x, y[, z] : array_like
            The vertex coordinate values as 1D arrays

        connect : ['strip'] | 'segments' or numpy.ndarray
            Determines which vertices are connected by lines
        
        color : str
            The color of the line

        symbol : str
            Marker symbol to use

        width : float
            The line width in px

        marker_size : float
           The marker size in px (if `0` no markers will be shown)
        """

        if z is None:
            self._configure(proj='2d')
            data = (x, y)
        else:
            self._configure(proj='3d')
            data = (x, y, z)

        line = scene.LinePlot(data=data, color=color, symbol=symbol,
                              width=width, marker_size=marker_size, connect=connect)

        self.visuals += [line]
        self.view.add(line)
        self.view.camera.set_range()

        return line

    def triplot(self, x, y, faces, color='black'):
        """
        Draw an unstructured 2-dimensional grid.

        Parameters
        ----------

        x, y, z : array_like
            The vertex coordinate values as 1D arrays
        
        faces : array_like
            The connectivity array of the triangular faces

        color : str
            The color to use
        """

        self._configure(proj='2d')

        return self.wireframe(x=x, y=y, z=np.zeros_like(x),
                              faces=faces, edge_color=color)
        
    def trisurface(self, x, y, z, faces, edges=None,
                   vertex_colors=None, face_colors=None, edge_color=None,
                   color=None, vmin=None, vmax=None):
        """
        Draw a surface consisting of triangular patches.

        Parameters
        ----------

        x, y, z : array_like
            The vertex coordinate values as 1D arrays
        
        faces : array_like
            The connectivity array of the triangular faces

        edges : array_like | None
            The connectivity array of the edges

        vertex_colors : array_like | None
            Colors to use for each vertex

        face_colors : array_like | None
            Colors to use for each face

        edge_color : str | tuple | None
            Color to use for the edges

        color : str | tuple | None
            Color to use

        vmin, vmax : float
            Min/Max values for the colormap
        """

        self._configure(proj='3d')
        
        vertices = np.column_stack([x, y, z]).astype(np.float32)
        faces = faces.astype(np.uint32)

        vmin = vmin if vmin is not None else z.min()
        vmax = vmax if vmax is not None else z.max()

        if all(c is None for c in [color, vertex_colors, face_colors]):
            if not np.allclose(vmin, vmax):
                t = (z - vmin) / (vmax - vmin)
            else:
                t = (z - 0.5*vmin)

            vertex_colors = self.cmap.map( np.expand_dims(t, axis=1) )

        mesh = visuals.SurfaceMesh(vertices=vertices, faces=faces, edges=edges,
                                   vertex_colors=vertex_colors, face_colors=face_colors,
                                   edge_color=edge_color, color=color)

        self.visuals += [mesh]
        self._autoscale()
        self.view.add(mesh)
        
        return mesh

    def wireframe(self, x, y, z, faces, edges=None, edge_color='black'):
        """
        Draw a wireframe outline of a mesh.

        Parameters
        ----------

        x, y, z : array_like
            The vertex coordinate values as 1D arrays
        
        faces : array_like
            The connectivity array of the triangular faces

        edges : array_like | None
            The connectivity array of the edges

        edge_color : str | tuple | None
            Color to use for the edges
        """

        self._configure(proj='3d')

        vertices = np.column_stack([x, y, z]).astype(np.float32)
        faces = faces.astype(np.uint32)

        mesh = visuals.WireframeMesh(vertices=vertices, faces=faces, edges=edges,
                                     edge_color=edge_color)

        self.visuals += [mesh]
        self._autoscale()
        self.view.add(mesh)

        return mesh

    def scatter2d(self, x, y, c=None, size=1, vertex_colors=None):
        """
        Plot scattered 2d data points.

        Parameters
        ----------

        x, y : array_like
            The vertex coordinate values as 1D arrays

        c : array_like
            Scalar values to use if vertex colors are not given

        size : float
            The point size in px

        vertex_colors : array_like | None
            Colors to use for each vertex

        color : tuple | None
            Color to use
        """
        self._configure(proj='2d')
        z = np.zeros_like(x)

        return self.scatter3(x=x, y=y, z=z, c=c, size=size,
                             vertex_colors=vertex_colors)

    def scatter3d(self, x, y, z, c=None, size=1,
                  colors=None, alpha=1.):
        """
        Plot scattered 3d data points.

        Parameters
        ----------

        x, y, z : array_like
            The vertex coordinate values as 1D arrays

        c : array_like
            Scalar values to use if `colors` are not given

        size : float
            The point size in px

        colors : str | tuple | array_like | None
            Colors used to draw each data point

        alpha : float
            The alpha value for the colors (opacity)
        """
        self._configure(proj='3d')

        pos = np.column_stack([x, y, z]).astype(np.float32)

        if colors is None:
            if c is not None:
                colors = self._get_vertex_colors(values=c, alpha=alpha)
            else:
                colors = np.repeat([(0,0,0,alpha)], np.size(vertices, axis=0), axis=0)
                
        scatter = visuals.ScatterPlot(pos=pos, size=size, colors=colors)

        self.visuals += [scatter]
        self._autoscale()
        self.view.add(scatter)

        return scatter
