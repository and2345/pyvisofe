"""
Provides some visuals that can be drawn in a scene.
"""

# IMPORTS
import numpy as np

from vispy import gloo, scene, visuals, io
from vispy.color import get_colormap

from .  import utils

class WireframeMeshVisual(visuals.MeshVisual):
    """
    Visual that displays a wireframe mesh.

    Parameters
    ----------

    vertices : array_like
        The vertex coordinates
    
    faces : array_like
        The connectivity array of the triangular faces

    edges : array_like | None
        The connectivity array of the edges

    edge_color : str | tuple | None
        Color to use for the edges
    """

    def __init__(self, vertices=None, faces=None, edges=None, edge_color='black'):

        # create wireframe mesh visual
        if edges is None and faces is not None:
            edges = utils.unique_rows(np.vstack(faces.take([[0, 1],
                                                            [0, 2],
                                                            [1, 2]], axis=1)))

        # initialize mesh visual
        visuals.MeshVisual.__init__(self, vertices=vertices, faces=edges.ravel(order='C'),
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
    """

    def __init__(self, vertices=None, faces=None, edges=None,
                 vertex_colors=None, face_colors=None, edge_color=None, color=None):

        if all(c is None for c in (vertex_colors, face_colors, color)):
            color = (0.5, 0.5, 1.0, 1.0)
        
        # create mesh and outline visuals
        self._mesh = visuals.MeshVisual(vertices=vertices, faces=faces,
                                        vertex_colors=vertex_colors,
                                        face_colors=face_colors,
                                        color=color,
                                        shading=None, mode='triangles')

        if edge_color is not None:
            self._outline = WireframeMesh(vertices=vertices, faces=faces,
                                          edges=edges, edge_color=edge_color)
        else:
            self._outline = visuals.MeshVisual()

        # initialize compound
        visuals.CompoundVisual.__init__(self, subvisuals=[self._mesh, self._outline])
            
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

class ScatterPlotVisual(visuals.Visual):
    """
    Visual the plot scattered data points.

    Parameters
    ----------

    vertices : array_like
        The vertex coordinates

    size : float
        The point size
    
    vertex_colors : array_like | None
        Colors to use for each vertex
    
    color : str | tuple | None
        Color to use
    """

    vertex_shader = """
varying vec4 v_color;

void main() {
    gl_Position = $transform(vec4($position, 1));
    gl_PointSize = $size;

    v_color = $color;
}
"""
    fragment_shader = """
varying vec4 v_color;

void main() {
    gl_FragColor = v_color;
}
"""
    """
    //gl_FragColor = $color;

    //vec2 pos = mod(gl_FragCoord.xy, vec2(50.0)) - vec2(25.0);
    //float sqdist = dot(pos, pos);
    //gl_FragColor = mix(vec4())

    //float d = 1 - length(gl_PointCoord - vec2(.5,.5)) / (sqrt(2)/2);
    //gl_FragColor = d * $color;
    //gl_FragColor.a = d;
    """

    def __init__(self, vertices=None, size=1.,
                 vertex_colors=None, color=None):

        if vertices is not None:
            vertices = np.asarray(vertices, dtype=np.float32)
            if np.size(vertices, axis=1) == 2:
                vertices = np.column_stack([vertices,
                                            np.zeros_like(vertices[:,0])])

        if vertex_colors is not None:
            vertex_colors = np.asarray(vertex_colors, dtype=np.float32)
                
        visuals.Visual.__init__(self, self.vertex_shader, self.fragment_shader)

        self.vbo = gloo.VertexBuffer(vertices)
        self.cbo = gloo.VertexBuffer(vertex_colors)

        self.shared_program.vert['position'] = self.vbo
        self.shared_program.vert['size'] = size
        self.shared_program.vert['color'] = self.cbo
        #self.shared_program.frag['color'] = vertex_colors

        self._bounds = None
        if vertices is not None:
            self._bounds = zip(vertices.min(axis=0), vertices.max(axis=0))

        self.set_gl_state(depth_test=True, blend=True,
                          blend_func=('src_alpha', 'one_minus_src_alpha'))
            
        self._draw_mode = 'points'

    def _compute_bounds(self, axis, view):
        if self._bounds is None:
            return None
        else:
            return self._bounds[axis]
        
    def _prepare_transforms(self, view):
        view.view_program.vert['transform'] = view.get_transform()
        
# these are the actual visuals that can be added to a scene
WireframeMesh = scene.visuals.create_visual_node(WireframeMeshVisual)
SurfaceMesh = scene.visuals.create_visual_node(SurfaceMeshVisual)
ScatterPlot = scene.visuals.create_visual_node(ScatterPlotVisual)
