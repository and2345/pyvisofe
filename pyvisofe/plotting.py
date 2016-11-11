# IMPORTS
from .canvas import Figure

# some wrappers for functional visualization
def plot(x, y, **kwargs):
    """
    See :py:func:`ViewWidget.plot`.
    """
    fig = Figure(show=True)
    ax = fig[0, 0]

    return ax.plot(x, y, **kwargs)

def triplot(x, y, triangles, **kwargs):
    """
    See :py:func:`ViewWidget.triplot`.
    """
    fig = Figure(show=True)
    ax = fig[0, 0]

    return ax.triplot(x, y, triangles, **kwargs)

def trisurface(x, y, z, triangles, **kwargs):
    """
    See :py:func:`ViewWidget.trisurface`.
    """
    fig = Figure(show=True)
    ax = fig[0, 0]

    return ax.trisurface(x, y, z, triangles, **kwargs)

def wireframe(x, y, z, triangles, **kwargs):
    """
    See :py:func:`ViewWidget.wireframe`.
    """
    fig = Figure(show=True)
    ax = fig[0, 0]

    return ax.wireframe(x, y, z, triangles, **kwargs)

# def heatmap(x, y, z, triangles, **kwargs):
#     """
#     See :py:func:`ViewWidget.heatmap`.
#     """

#     # instanciate canvas and get ViewWidget axes
#     fig = Figure(show=True)
#     ax = fig[0, 0]

#     # use ViewWidget's method
#     return ax.heatmap(x=x, y=y, z=z, triangles=triangles, **kwargs)

