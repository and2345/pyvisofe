"""
Provides a canvas to draw on.
"""

# IMPORTS
import numpy as np

from vispy import scene
from vispy.color import get_colormap

from . import widget

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
        self._grid._default_class = widget.ViewWidget

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
            
