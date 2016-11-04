"""
Provides some visuals that can be drawn in a scene.
"""

# IMPORTS
import numpy as np
from vispy import gloo, scene, visuals, io
from vispy.color import get_colormap

from ..utils import unique_rows

class WireframeMeshVisual(visuals.MeshVisual):
    """
    Visual that displays a wireframe mesh.

    Parameters
    ----------

    vertices : array_like
        The vertex coordinates
    
    faces : array_like
        The face connectivity

    edges : array_like
        The edge connectivity

    edge_color : str | vispy.color.Color
        The color to use for the faces outline
    """

    def __init__(self, vertices=None, faces=None, edges=None, edge_color='black', **kwargs):

        # create wireframe mesh visual
        if edges is None and faces is not None:
            edges = unique_rows(np.vstack(faces.take([[0, 1],
                                                      [0, 2],
                                                      [1, 2]], axis=1)))

        # initialize mesh visual
        visuals.MeshVisual.__init__(self, vertices=vertices, faces=edges.ravel(),
                                    color=edge_color, mode='lines')

    @property
    def edge_color(self):
        """
        The color of the edges.
        """
        return self.color

    @edge_color.setter
    def edge_color(self, color):
        self.color = color

class SurfaceMeshVisual(visuals.CompoundVisual):
    """
    Visual that displays a surface mesh.

    Parameters
    ----------

    vertices : array_like
        The vertex coordinates
    
    faces : array_like
        The face connectivity

    vertex_colors : array_like
        Colors to use for each vertex
    
    face_colors : array_like
        Colors to use for each face

    color : str | vispy.color.Color
        The color to use for the faces

    edge_color : str | vispy.color.Color
        The color to use for the faces outline

    edges : array_like
        The edge connectivity
    """

    def __init__(self, vertices=None, faces=None, edges=None,
                 vertex_colors=None, face_colors=None, edge_color=None, color=None,
                 **kwargs):

        if vertex_colors is None and face_colors is None and color is None:
            color = (0.5, 0.5, 1.0, 1.0)
        
        # create mesh and outline visuals
        self._mesh = visuals.MeshVisual(vertices=vertices, faces=faces,
                                        vertex_colors=vertex_colors,
                                        face_colors=face_colors,
                                        color=color)

        if edge_color is not None:
            self._outline = WireframeMesh(vertices=vertices, faces=faces,
                                               edge_color=edge_color, edges=edges)
        else:
            self._outline = visuals.MeshVisual()

        # initialize compound
        visuals.CompoundVisual.__init__(self, subvisuals=[self._mesh, self._outline], **kwargs)
            
        # set polygon offset to make outlines visible
        self._mesh.set_gl_state(polygon_offset_fill=True,
                                polygon_offset=(1, 1),
                                depth_test=True)
    
    @property
    def mesh(self):
        """
        The vispy.visuals.MeshVisual that used to draw the filled triangles.
        """
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        self._mesh = mesh

    @property
    def mesh_data(self):
        """
        The mesh data.
        """
        return self._mesh.mesh_data

    def mesh_data_changed(self):
        return self._mesh.mesh_data_changed()

    @property
    def outline(self):
        """
        The vispy.visuals.MeshVisual that used to draw the triangle outlines.
        """
        return self._outline

    @outline.setter
    def outline(self, outline):
        self._outline = outline

# these are the actual visuals that can be added to a scene
WireframeMesh = scene.visuals.create_visual_node(WireframeMeshVisual)
SurfaceMesh = scene.visuals.create_visual_node(SurfaceMeshVisual)
