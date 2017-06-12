# !/usr/bin/env
# Winged edge polyhedron model


import math
from random import randrange
import copy
import warnings
import itertools

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# import scipy as sc
import mpl_toolkits.mplot3d as a3
from matplotlib.patches import Polygon as mplPoly
from matplotlib.collections import PatchCollection


COMPARISON_EPSILON = 0.0001


class vertex:
    """Class representing a vertex in a three-dimensional space"""

    def __init__(self, x=0, y=0, z=0):
        """ float * float * float / (float * float * float) -> vertex
        Builds a vertex using its coordinates.
        If three arguments are given, they are interpreted as the x, y and z
            coordinates of the vertex.
        If one argument is given, it is interpreted as a tuple of the three
            coordinates."""
        if type(x) is tuple:
            self.x = x[0]
            self.y = x[1]
            self.z = x[2]
        elif type(x) is vector:
            self.x = x.x
            self.y = x.y
            self.z = x.z
        else:
            self.x = x
            self.y = y
            self.z = z

    def __str__(self):
        return ("({}, {}, {})".format(self.x, self.y, self.z) +
                ("" if not self.isInterior() else " interior"))

    def coords(self):
        """ vertex -> float**3
        Returns the tuple of coordinates of <self>."""
        return (self.x, self.y, self.z)

    def markInt(self):
        """ vertex -> None
        Marks the vertex as interior."""
        self.interior = True
        # print('int')

    def isInterior(self):
        """ vertex -> bool"""
        try:
            return self.interior
        except AttributeError:
            return False

    def dist(self, other):
        """ vertex * vertex -> float
        Returns the distance between the two vertices."""
        return math.sqrt((self.x - other.x)**2 +
                         (self.y - other.y)**2 +
                         (self.z - other.z)**2)

    def copy(self):
        """ vertex -> vertex
        Returns a copy of <self>."""
        return vertex(self.x, self.y, self.z)

    def markFTraced(self, f):
        """ vertex * face -> None
        Marks <self> as having been traced for the face <f>."""
        try:
            self.fTraced.add(f)
        except AttributeError:
            self.fTraced = set()
            self.fTraced.add(f)

    def isMarkedFor(self, f):
        try:
            return f in self.fTraced
        except AttributeError:
            return False


class edge:
    """Class representing an oriented edge in a three-dimensional space"""

    def __init__(self, nVertex, pVertex=None, pFace=None, nFace=None):
        """ vertex * vertex * face * face -> edge
        Builds an oriented edge using its adjacent vertices and faces.
        It is okay for the faces to be None, if the edge is not linked.
        However None vertices are not justifiable and shall not be used."""
        if type(nVertex) is vector:
            self.pvt =  vertex(nVertex)
            self.nvt = vertex(0, 0, 0)
        else:
            self.pvt = pVertex  # positive vertex
            self.nvt = nVertex  # negative vertex
        self.pFace = pFace  # positive face
        self.nFace = nFace  # negative face

    def linkFace(self, f):
        """ edge * face -> ()
        Links the face <f> to the edge <self>."""
        if self.pFace is None:
            self.pFace = f
        elif self.nFace is None:
            self.nFace = f
        else:
            raise ValueError('Edge is already linked to two faces')

    def __str__(self):
        return "{} -> {}".format(str(self.nvt), str(self.pvt))

    def invert(self):
        """ edge -> ()
        Reverses the direction of the edge (exchanges the pvt with nvt
        and the pFace with the nFace)."""
        tmp = self.pvt
        self.pvt = self.nvt
        self.nvt = tmp
        tmp = self.pFace
        self.pFace = self.nFace
        self.nFace = tmp

    def faceIntersection(self, face):
        """ edge * face -> float**3/bool
        Returns the vertex at the intersection of <self> and <face>
            (if it exists).
        If they do not intersect, returns False."""
        # Test wether both ends of the edge are in the same half space
        # (relative to <face>'s plane).
        normal = face.normalVect()
        v0 = vector(face.vertices[0])
        vp = vector(self.pvt)
        vn = vector(self.nvt)
        p = normal.dotProduct(vp - v0) * normal.dotProduct(vn - v0)
        if p > 0:
            return False
        elif abs(p) <= COMPARISON_EPSILON or abs(normal.dotProduct(vp - vn) / (normal.norm() * (vp - vn).norm())) <= COMPARISON_EPSILON:
            # print('ah')
            return False
        else:
            interVect = vn + (normal.dotProduct(v0 - vn) /
                              normal.dotProduct(vp - vn)) * (vp - vn)
            lastCross = ((vector(face.vertices[-1]) - interVect) *
                         (vector(face.vertices[0]) - interVect))
            for i in range(len(face.vertices)):
                cross = ((vector(face.vertices[i]) - interVect) *
                         (vector(face.vertices[(i + 1) % len(face.vertices)]) -
                          interVect))
                p = cross.dotProduct(lastCross)
                if p < 0:
                    return False
                elif p == 0 and cross.norm() != 0:
                    if cross.norm() > COMPARISON_EPSILON:
                        warnings.warn("Cross product's norm is very low")
                    lastCross = cross
            return interVect.coords()

    def polyhedronIntersection(self, poly):
        """ edge * polyhedron -> vertex list
        Returns the list of vertices at the intersections of
            <self> and <poly>."""
        return list(filter(lambda x: x is not False,
                           (self.faceIntersection(f) for f in poly.faces)))

    def other(self, fv):
        """ edge * face/vertex -> face/vertex
        Returns the face/vertex opposite from <fv> with respect to <self>"""
        if self.pvt == fv:
            return self.nvt
        elif self.nvt == fv:
            return self.pvt
        elif self.pFace == fv:
            return self.nFace
        elif self.nFace == fv:
            return self.pFace

    def markIntersectedWith(self, face):
        """ edge * face -> None
        Marks <self> as having already been intersected with <face>."""
        try:
            self.hasIntersected.add(face)
        except AttributeError:
            self.hasIntersected = set()
            self.hasIntersected.add(face)

    def hasIntersectedWith(self, f):
        """ edge * face -> bool
        Returns wether <self> has already been intersected with <f>."""
        try:
            return f in self.hasIntersected
        except AttributeError:
            return False


class face:
    """Class representing an oriented face in a three-dimensional space"""

    def __init__(self, vertices):
        """ vertex list -> face
        Builds an oriented face from a list of its vertices.
        Coplanarity of the vertices will not be checked."""
        self.vertices = vertices

    def __str__(self):
        return "({})".format(' - '.join(str(v) for v in self.vertices))

    def invert(self):
        """ face -> ()
        Flips the face by reversing the order of the vertices
            in the vertices list."""
        self.vertices.reverse()

    def normalVect(self, n=2):
        """ face * int -> vector
        Returns the vector normal to the face (normalized)."""
        L = len(self.vertices)
        normals = []
        while len(normals) < n:
            j = randrange(L)
            v0 = vector(self.vertices[j].coords())
            v1 = vector(self.vertices[int(j + L / 3) % L].coords())
            v2 = vector(self.vertices[int(j + 2 * L / 3) % L].coords())
            try:
                normals.append(((v1 - v0) * (v2 - v0)).normalize())
            except ValueError:
                pass
        return (1 / len(normals)) * sum(normals, vector(0, 0, 0))

    def copy(self):
        """ face -> face
        Returns a deep copy of <self>."""
        newVertices = [v.copy() for v in self.vertices]
        return face(newVertices)

    def addVertexBetween(self, v1, nV, v2):
        """ face * vertex * vertex * vertex -> None
        Adds a new vertex in <self>.vertices between the two vertices given."""
        n = len(self.vertices)
        if nV not in self.vertices:
            for i in range(n):
                if self.vertices[i] == v1:
                    if self.vertices[(i + 1) % n] == v2:
                        self.vertices.insert((i + 1) % n, nV)
                        return nV
                    elif self.vertices[(i - 1) % n] == v2:
                        self.vertices.insert(i, nV)
                        return nV
                    else:
                        # print(v1, v2, self)
                        raise ValueError('Vertices given are not '
                                         'adjacent in the face')
            # # print(v1, self)
            raise ValueError('<v1> is not even in the face')

    def isOnInteriorSide(self, v):
        """ face * vertex -> bool
        Returns wether <v> is on the interior side of <self>."""
        n = self.normalVect()
        return n.dotProduct(vector(self.vertices[0]) - vector(v)) > 0

    def edges(self):
        """ face -> edges list
        Returns the list of oriented edges constituting the face."""
        return [edge(self.vertices[i - 1], self.vertices[(i)]) for i in range(-1, len(self.vertices))]

    def containsEdge(self, e):
        """ face * edge -> bool
        Returns wether <self> contains <e> in its list of vertices."""
        return any(e.nvt in [self.vertices[i-2], self.vertices[i]] and self.vertices[i-1] == e.pvt for i in range(len(self.vertices)))

    def hasNoDoubleVertices(self):
        """ face -> bool
        Ensures <self> has no double vertices."""
        assert all(self.vertices.count(v) == 1 for v in self.vertices)
        return (all(all(v1 == v2 or v1.dist(v2) > COMPARISON_EPSILON for v2 in self.vertices)
                   for v1 in self.vertices) and
                all(self.vertices.count(v) == 1 for v in self.vertices))


class polyhedron:
    """Class representing a polyhedron in a three-dimensional space"""

    def __init__(self, vertices, edges, faces):
        """vertex list * edge list * face list -> polyhedron
        Builds a polyhedron using lists of its vertices, edges and faces.
        The coherence of the polyhedron is not verified."""
        self.vertices = vertices
        self.edges = edges
        self.faces = faces

    def __str__(self):
        return """polyhedron( \n\tVertices: {0}
                              \n\tEdges: {1}
                              \n\tFaces: {2}\n)"""\
                .format('\n\t\t'.join(str(v) for v in self.vertices),
                        '\n\t\t'.join(str(e) for e in self.edges),
                        '\n\t\t'.join(str(f) for f in self.faces))

    def makeBasicPolyhedron():
        """ None -> polyhedron
        Returns a basic polyhedron with one vertex (at (0,0,0))
            and one face with the vertex.
        This is the most basic polyhedron which satisfies Euler's
            polyhedral formula.
        Corresponds to Baumgart's MKBFV function."""
        v = vertex(0, 0, 0)
        return polyhedron([v], [], [face([v])])

    def getVertex(self, x, y, z, epsilon=COMPARISON_EPSILON):
        """ polyhedron * float * float * float * float
        Returns (if  it exists) the vertex at <x>, <y> and <z>
            coordinates, give or take epsilon."""
        for v in self.vertices:
            if (v.x - x)**2 + (v.y - y)**2 + (v.z - z)**2 <= epsilon**2:
                return v
        raise ValueError('No vertex found')

    def getEdge(self, v1, v2):
        """ polyhedron * vertex * vertex -> edge
        Returns (if it exists) the edge between <v1> and <v2>."""
        for e in self.edges:
            if (e.pvt, e.nvt) in [(v1, v2), (v2, v1)]:
                return e
        raise ValueError('No edge found')

    def getFace(self, vertices):
        """ polyhedron * vertex list -> face
        Returns (if  it exists) the face that has <vertices> as vertices."""
        for f in self.faces:
            if f.vertices == vertices:
                return f
        raise ValueError('No face found')

    def containsEdge(self, v1, v2):
        """ polyhedron * vertex * vertex -> bool
        Returns wether or not the edge between <v1> and <v2> exists in <self>.
        NOTE: this checks for identity between vertices, NOT for equality.
            If an edge exists in <self> at the exact coordinates
            of <v1> and <v2> but its end vertices are not the same object
            as <v1> and <v2>, the output will be False"""
        for e in self.edges:
            if (e.pvt, e.nvt) in [(v1, v2), (v2, v1)]:
                return True
        return False

    def addFace(self, vertices, bypassCheck=False):
        """ polyhedron * vertex list -> face
        Creates a face and adds it to the polyhedron while
            keeping it coherent and valid, and returns the face.
        Missing vertices will be created and added.
        If the face cannot be added (e.g. if one of the
            edges of the vertices given already has two faces linked),
            a ValueError exception will be raised."""
        try:
            if bypassCheck:
                raise ValueError
            return self.getFace(vertices)
        except ValueError:
            if any(vertices.count(v) > 1 for v in vertices):
                raise ValueError('The face given is invalid: '
                                 'two or more vertices are identical')
            newF = face(vertices)
            self.faces.append(newF)
            for i in range(len(vertices)):
                try:
                    e = self.getEdge(vertices[i],
                                     vertices[(i + 1) % len(vertices)])
                except ValueError:
                    e = self.addEdge(vertices[i],
                                     vertices[(i + 1) % len(vertices)])
                e.linkFace(newF)
            return newF

    def addEdge(self, nVertex, pVertex=None):
        """ polyhedron * edge/vertex * vertex -> edge
        Creates an edge between <nVertex> and <pVertex> and adds it to <self>,
            if it did not exist previously in <self>, and returns it.
        If it has to be created,
            adjacents faces for the new edge are set to None."""
        if not self.containsEdge(nVertex, pVertex):
            if type(nVertex) is edge:
                self.edges.append(nVertex)
                self.addVertex(nVertex.pvt)
                self.addVertex(nVertex.nvt)
                return nVertex
            else:
                self.addVertex(nVertex)
                self.addVertex(pVertex)
                newE = edge(nVertex, pVertex)
                self.edges.append(newE)
                return newE

        else:
            return self.getEdge(nVertex, pVertex)

    def addVertex(self, arg1, arg2=0, arg3=0):
        """ polyhedron * (float/vector/(float * float * float)) *
            float * float -> vertex
        Creates a vertex at given coordinates, adds it to <self>,
            and returns it."""
        if type(arg1) is tuple:
            x, y, z = arg1
        elif type(arg1) is vector:
            x, y, z = arg1.coords()
        elif type(arg1) is float or type(arg1) is int:
            x, y, z = arg1, arg2, arg3
        elif type(arg1) is vertex:
            try:
                newV = self.getVertex(arg1.x, arg1.y, arg1.z, COMPARISON_EPSILON)
            except ValueError:
                newV = arg1
                self.vertices.append(arg1)
            return newV
        else:
            raise ValueError('bad argument type: ' + str(type(arg1)))
        try:
            newV = self.getVertex(x, y, z, COMPARISON_EPSILON)
        except ValueError:
            newV = vertex(x, y, z)
            self.vertices.append(newV)
        return newV

    def nnw(self, e):
        """ polyhedron * edge -> edge
        Returns the negative-negative wing of the edge e
            (the one that borders the negative face of <e>,
            and touches the negative vertex of <e>).
        Equivalent to NCW in Baumgart's conventions."""
        vList = self.nFace.vertices
        n = len(vList)
        i = vList.index(self.nvt)
        if vList[(i - 1) % n] == self.pvt:
            return self.getEdge(self.nvt, vList[(i + 1) % n])
        else:
            return self.getEdge(self.nvt, vList[(i - 1) % n])

    def pnw(self, e):
        """ polyhedron * edge -> edge
        Returns the positive-negative wing of the edge e
            (the one that borders the positive face of <e>,
            and touches the negative vertex of <e>).
        Equivalent to PCCW in Baumgart's conventions."""
        vList = self.pFace.vertices
        n = len(vList)
        i = vList.index(self.nvt)
        if vList[(i - 1) % n] == self.pvt:
            return self.getEdge(self.nvt, vList[(i + 1) % n])
        else:
            return self.getEdge(self.nvt, vList[(i - 1) % n])

    def npw(self, e):
        """ polyhedron * edge -> edge
        Returns the negative-positive wing of the edge e
            (the one that borders the negative face of <e>,
            and touches the positive vertex of <e>).
        Equivalent to NCCW in Baumgart's conventions."""
        vList = self.nFace.vertices
        n = len(vList)
        i = vList.index(self.pvt)
        if vList[(i - 1) % n] == self.nvt:
            return self.getEdge(self.pvt, vList[(i + 1) % n])
        else:
            return self.getEdge(self.pvt, vList[(i - 1) % n])

    def ppw(self, e):
        """ polyhedron * edge -> edge
        Returns the positive-positive wing of the edge e
            (the one that borders the positive face of <e>,
            and touches the positive vertex of <e>).
        Equivalent to PCW in Baumgart's conventions."""
        vList = self.pFace.vertices
        n = len(vList)
        i = vList.index(self.pvt)
        if vList[(i - 1) % n] == self.nvt:
            return self.getEdge(self.pvt, vList[(i + 1) % n])
        else:
            return self.getEdge(self.pvt, vList[(i - 1) % n])

    def nextEdgeCW(self, e, fv):
        """ polyhedron * edge * face/vertex -> edge
        Returns the next edge bordering <fv>,
            rotating clockwise and starting from <e>.
        Equivalent to ECW in Baumgart's conventions"""
        if self.pFace == fv:
            return self.ppw(e)
        elif self.nFace == fv:
            return self.nnw(e)
        elif self.pvt == fv:
            return self.npw(e)
        elif self.nvt == fv:
            return self.pnw(e)
        else:
            raise ValueError('The edge given does not border '
                             'the face/vertex given')

    def nextEdgeCCW(self, e, fv):
        """ polyhedron * edge * face/vertex -> edge
        Returns the next edge bordering <fv>,
            rotating counter-clockwise and starting from <e>.
        Equivalent to ECCW in Baumgart's conventions"""
        if self.pFace == fv:
            return self.pnw(e)
        elif self.nFace == fv:
            return self.npw(e)
        elif self.pvt == fv:
            return self.ppw(e)
        elif self.nvt == fv:
            return self.nnw(e)
        else:
            raise ValueError('The edge given does not border '
                             'the face/vertex given')

    def nextFaceCW(self, e, v):
        """ polyhedron * edge * vertex -> face
        Returns the next face bordering <v>,
            rotating clockwise with respect to <v> and starting from <e>"""
        if self.pvt == v:
            return self.nFace
        elif self.nvt == v:
            return self.pFace
        else:
            raise ValueError('The vertex given does not border the edge given')

    def nextFaceCCW(self, e, v):
        """ polyhedron * edge * vertex -> face
        Returns the next face bordering <e> and <v>,
            rotating counter-clockwise with respect to <v>
            and starting from <e>"""
        if self.pvt == v:
            return self.pFace
        elif self.nvt == v:
            return self.nFace
        else:
            raise ValueError('The vertex given does not border the edge given')

    def nextVertexCW(self, e, f):
        """ polyhedron * edge * face -> vertex
        Returns the next vertex bordering <e> and <f>,
            rotating clockwise with respect to <f> and starting from <e>"""
        if self.pFace == f:
            return self.pvt
        elif self.nFace == f:
            return self.nvt
        else:
            raise ValueError('The face given does not border the edge given')

    def nextVertexCCW(self, e, f):
        """ polyhedron * edge * face -> vertex
        Returns the next vertex bordering <e> and <f>,
            rotating counter-clockwise with respect to <f>
            and starting from <e>"""
        if self.pFace == f:
            return self.nvt
        elif self.nFace == f:
            return self.pvt
        else:
            raise ValueError('The face given does not border the edge given')

    def linked(self, el1, el2):
        """ polyhedron * vertex/edge/face * vertex/edge/face -> vertex/edge/bool
        Returns wether to <el1> and <el2> are linked or not and, in some cases,
            the object that links them."""
        if type(el1) is face and type(el2) is face:
            n1 = len(el1.vertices)
            n2 = len(el2.vertices)
            for i in range(n1):
                for j in range(n2):
                    if el1.vertices[i] == el2.vertices[j]:
                        if (el1.vertices[(i + 1) % n1]
                                in [el2.vertices[(j + 1) % n2],
                                    el2.vertices[(j - 1) % n2]]):
                            return self.getEdge(el1.vertices[i % n1],
                                                el1.vertices[(i + 1) % n1])
                        elif (el1.vertices[(i - 1) % n1]
                                in [el2.vertices[(j + 1) % n2],
                                    el2.vertices[(j - 1) % n2]]):
                            return self.getEdge(el1.vertices[i % n1],
                                                el1.vertices[(i - 1) % n1])
            return False

        elif type(el1) is face and type(el2) is edge:
            return el1 in [el2.pFace, el2.nFace]
        elif type(el1) is edge and type(el2) is face:
            return self.linked(el2, el1)

        elif type(el1) is edge and type(el2) is edge:
            if el1.pvt in [el2.pvt, el2.nvt]:
                return el1.pvt
            elif el1.nvt in [el2.pvt, el2.nvt]:
                return el1.nvt
            else:
                return False

        elif type(el1) is edge and type(el2) is vertex:
            return el2 in [el1.pvt, el1.nvt]
        elif type(el1) is vertex and type(el2) is edge:
            return self.linked(el2, el1)

        elif type(el1) is vertex and type(el2) is vertex:
            try:
                return self.getEdge(el1, el2)
            except ValueError:
                return False

        else:
            raise TypeError('Invalid arguments types')

    def evert(self):
        """ polyhedron -> ()
        Turns <self> inside out (inverts all of its edges and faces)."""
        for e in self.edges:
            self.invert()
        for f in self.faces:
            f.invert()

    def edgeSplit(self, e, newCoords=None):
        """ polyhedron * edge [* float**3] -> vertex * edge
        Splits <e> into two edges, one from <e>.nvt to <newVertex>,
            the other from <newVertex> to <e>.pvt.
        Returns the middle vertex and the new edgself.
        Default for <newVertex> is the middle point between <e>.nvt <e>.pvt.
        Corresponds to Baumgart's eulerian primitive ESPLIT."""
        assert (e.pFace is None or e.pFace.containsEdge(e)) and (e.nFace is None or e.nFace.containsEdge(e))
        assert all(all(v in self.vertices for v in f.vertices) for f in self.faces)
        if newCoords is None:
            newVertex = self.addVertex((e.nvt.x + e.pvt.x) / 2,
                                       (e.nvt.y + e.pvt.y) / 2,
                                       (e.nvt.z + e.pvt.z) / 2)
            assert newVertex in self.vertices
        else:
            newVertex = self.addVertex(newCoords)
            assert newVertex in self.vertices
        assert newVertex in self.vertices
        assert all(all(v in self.vertices for v in f.vertices) for f in self.faces)
        newEdge = self.addEdge(newVertex, e.pvt)
        assert newVertex in self.vertices
        assert all(all(v in self.vertices for v in f.vertices) for f in self.faces)
        e.pFace.addVertexBetween(e.nvt, newVertex, e.pvt)
        assert newVertex in self.vertices
        assert all(all(v in self.vertices for v in f.vertices) for f in self.faces)
        e.nFace.addVertexBetween(e.nvt, newVertex, newEdge.pvt)
        assert all(all(v in self.vertices for v in f.vertices) for f in self.faces)
        e.pvt = newVertex
        assert all(all(v in self.vertices for v in f.vertices) for f in self.faces)
        newEdge.pFace = e.pFace
        assert all(all(v in self.vertices for v in f.vertices) for f in self.faces)
        newEdge.nFace = e.nFace
        assert all(all(v in self.vertices for v in f.vertices) for f in self.faces)
        return (newVertex, newEdge)

    def faceSplit(self, f, v1, v2):
        """ polyhedron * face * vertex * vertex -> edge * face
        Splits <f> into two faces and creates an edge from <v1>
            and <v2> to separate them.
        Returns the separating edge and the new face.
        Corresponds to Baumgart's eulerian primitive MKFself."""
        if self.containsEdge(v1, v2):
            raise ValueError('There is already an edge between'
                             'the vertices given.'
                             'Is your polyhedron really one ?')
        try:
            i1 = f.vertices.index(v1)
            i2 = f.vertices.index(v2)
        except ValueError:
            raise ValueError('Vertices given are not on the face given')
        f.vertices, newVertices = splitList(f.vertices, i2, i1)
        newEdge = self.addEdge(v2, v1)
        newFace = self.addFace(newVertices)
        newEdge.pFace = f
        newEdge.nFace = newFace
        return (newEdge, newFace)

    def makeEdgeVertex(self, f, v):
        """ polyhedron * face * vertex -> edge * vertex
        Creates a dangling wire at vertex <v> around face <f>.
        Returns the dangling edge and vertex.
        The dangling vertex is created at (0,0,0).
        Corresponds to baumgart's Eulerian primitive MKEV."""
        newV = self.addVertex(0, 0, 0)
        newE = self.addEdge(v, newV)
        newE.pFace = f
        return (newE, newV)

    def translateSweep(self, vect, epsilon=COMPARISON_EPSILON):
        """ polyhedron * vector -> ()
        Sweeps the polyhedron by translating it following <vect>:
            vertices are transformed to edges, edges to quadrilaterals,
            faces to prism.
        If a sweeped vertex is to be created at the location of an old one
            (give or take <epsilon>), it will not be created.
            Set <epsilon> to negative to force the creation of all vertices.
        Sweeping a multi-faced polyhedron is not advised.
        Corresponds to one case of Baumgart's routine SWEEP,
            the other case being covered by the rotateSweep method."""
        return self.sweep(lambda c: (c[0] + vect.x,
                                     c[1] + vect.y,
                                     c[2] + vect.z),
                          epsilon)

    def rotateSweep(self, rotVect, rotCenter, epsilon=COMPARISON_EPSILON):
        """ polyhedron * vector * (float * float * float) * float -> ()
        Sweeps the polyhedron by rotating it around <rotVect> and <rotCenter>:
            vertices are transformed to edges,
            edges to quadrilaterals, faces to prism.
        <rotVect> is a rotation vector, given in radian.
        <rotCenter> is a tuple of 3 coordinates.
        If a sweeped vertex is to be created at the location of an old one
            (give or take <epsilon>), it will not be created.
            Set <epsilon> to negative to force the creation of all vertices.
        Sweeping a multi-faced polyhedron is not advised.
        Corresponds to one case of Baumgart's routine SWEEP,
            the other case being covered by the translateSweep method."""
        rotn = (1 / rotVect.norm()) * rotVect

        def rot(coords):
            oldVect = vector(coords) - vector(rotCenter)
            newVect = (math.cos(rotVect.norm()) * oldVect +
                       math.sin(rotVect.norm()) * rotn * oldVect +
                       ((1 - math.cos(rotVect.norm())) *
                        (rotn.dotProduct(oldVect)) * rotn))
            return (rotCenter[0] + newVect.x,
                    rotCenter[1] + newVect.y,
                    rotCenter[2] + newVect.z)
        return self.sweep(rot, epsilon)

    def sweep(self, f, epsilon=COMPARISON_EPSILON):
        """ polyhedron * vector *
                ((float * float * float) -> (float * float * float)) *
                float
            -> ()
        Sweeps the polyhedron by applying <f> to its points:
            vertices are transformed to edges,
            edges to quadrilaterals, faces to prism.
        If a sweeped vertex is to be created at the location of an old one
            (give or take <epsilon>), it will not be created.
            Set <epsilon> to negative to force the creation of all vertices.
        Sweeping a multi-faced polyhedron is not advised.
        Corresponds loosely to Baumgart's routine SWEEP."""
        nv = len(self.vertices)
        newVertices = []
        oldEdges = self.edges.copy()
        oldFaces = self.faces.copy()
        excludeV = []
        for i in range(nv):
            v = self.vertices[i]
            newVertex = vertex(f(v.coords()))
            try:
                collidingV = self.getVertex(newVertex.x,
                                            newVertex.y,
                                            newVertex.z,
                                            epsilon)
                newVertices.append(collidingV)
                excludeV.append(collidingV)
            except ValueError:
                newVertices.append(newVertex)
                # newVertices[i] now corresponds to self.vertices[i] sweeped
            self.addEdge(v, newVertices[i])
        for e in oldEdges:
            v1 = e.nvt
            v2 = e.pvt
            tv1 = newVertices[self.vertices.index(v1)]
            tv2 = newVertices[self.vertices.index(v2)]
            self.addEdge(tv1, tv2)
            try:
                self.addFace([v1, v2, tv2, tv1])
            except ValueError:
                pass
        for f in oldFaces:
            try:
                self.addFace([newVertices[self.vertices.index(v)]
                              for v in f.vertices])
            except ValueError:
                pass
        self.vertices.extend(filter(lambda x: x not in excludeV, newVertices))

    def plot(self, plotEdges=False, emphaseEdges=[], col=('b', 'k', 'r'), lims=None, ort=False):
        """ polyhedron -> ()
        Plots the polyhedron using matplotlib."""
        ax = a3.Axes3D(plt.figure())
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.dist = 30
        ax.azim = -140
        if lims is None:
            lims = [0, 0, 0]
            lims[0] = [min(v.x for v in self.vertices),
                       max(v.x for v in self.vertices)]
            lims[1] = [min(v.y for v in self.vertices),
                       max(v.y for v in self.vertices)]
            lims[2] = [min(v.z for v in self.vertices),
                       max(v.z for v in self.vertices)]
        if ort:
            ma = max(lims[i][1] for i in range(3))
            mi = min(lims[i][0] for i in range(3))
            lims = [[mi, ma]] * 3
        ax.set_xlim(lims[0])
        ax.set_ylim(lims[1])
        ax.set_zlim(lims[2])
        for f in self.faces:
            face = a3.art3d.Poly3DCollection([[v.coords()
                                               for v in f.vertices]])
            ax.add_collection3d(face)
            face.set_facecolor(col[0])
            face.set_edgecolor(col[1])
        if plotEdges or len(emphaseEdges)>0:
            for e in self.edges:
                edge = a3.art3d.Poly3DCollection([[e.nvt.coords(),
                                                   e.pvt.coords()]])
                ax.add_collection3d(edge)
                if e in emphaseEdges:
                    edge.set_edgecolor(col[2])
                else:
                    edge.set_edgecolor(col[1])
        plt.show()

    def apply(self, f):
        """ polyhedron * ((float * float * float) -> (float * float * float)) -> ()
        Applies <f> to all of <self>'s vertices."""
        for v in self.vertices:
            v.x, v.y, v.z = f(v.coords())

    def dilate(self, coef):
        """ polyhedron * float -> polyhedron
        Multiplies the coordinates of the vertices by coef."""
        self.apply(lambda c: (coef * vector(c)).coords())

    def translate(self, vect):
        """ polyhedron * vector -> ()
        Translates the <self> following <vect>."""
        self.apply(lambda c: (vector(c) + vect).coords())

    def union(self, other):
        """ polyhedron * polyhedron -> ()
        Merges <self> and <other> into one polyhedron."""
        self.vertices.extend(other.vertices)
        self.edges.extend(other.edges)
        self.faces.extend(other.faces)
        return self

    def intersection(self, other, bypassInclusionTest=False):
        """polyhedron * polyhedron -> polyhedron
        Computes and returns the polyhedron of the
            intersection of <self> and <other>."""
        # poly1 = copy.deepcopy(self)
        # poly2 = copy.deepcopy(other)
        poly1 = self
        poly2 = other
        noneF = []
        assert self.pnFacesInPoly()
        polyInter = polyIntersection.fromPolyhedron(poly1, poly2)
        assert all(all(v in polyInter.poly1.vertices for v in f.vertices) for f in polyInter.poly1.faces)
        assert all(all(v in polyInter.poly2.vertices for v in f.vertices) for f in polyInter.poly2.faces)
        # This handles the case of one polyhedron inside the other
        if not bypassInclusionTest:
            if len(polyInter.inter1) + len(polyInter.inter2) == 0:
                if all(f.isOnInteriorSide(poly1.vertices[0]) for f in poly2.faces):
                    return poly1
                elif all(f.isOnInteriorSide(poly2.vertices[0]) for f in poly1.faces):
                    return poly2
                else:
                    return polyhedron([], [], [])
        assert self.pnFacesInPoly() and self.pnFacesInPoly()
        # print('step 1')
        assert polyInter.nonDoubleIntersectors()
        polyInter.createIntersectorVertices()
        # print('step 2')
        polyInter.buildEdges()
        # print('step 3')
        polyInter.buildFaces()
        return polyInter.result

    def silhouettePhoto(self, focalPoint, focalDist, angles):
        """ polyhedron * float**3 * float * float**3 -> silhouette
        Returns the silhouette of <self>, seen from <focalPoint> with a
            focal of <focalDist>, and angles of shot specified by <angle>"""
        sil = silhouette(focalPoint, focalDist, angles, [])
        sil.edgesTaken = []
        X, Y, Z = baseFromAngles(angles)
        fVect = vector(focalPoint)
        fVert = vertex(focalPoint)
        for e in self.edges:
            vertOfP = next(filter(lambda x: x not in [e.pvt, e.nvt],
                                  e.pFace.vertices))
            vertOfN = next(filter(lambda x: x not in [e.pvt, e.nvt],
                                  e.nFace.vertices))
            equ = (X.x, X.y, X.z, focalDist + fVect.dotProduct(X))
            projOfP = vector(planeLineIntersect(focalPoint, vertOfP, equ))
            projOfN = vector(planeLineIntersect(focalPoint, vertOfN, equ))
            pvtProj = vector(planeLineIntersect(focalPoint, e.pvt, equ))
            nvtProj = vector(planeLineIntersect(focalPoint, e.nvt, equ))
            # # print((projOfP - fVect), vector(vertOfP) - fVect)
            # assert ((projOfP - fVect) * vector(vertOfP) - fVect).norm() < COMPARISON_EPSILON
            nvtReduced = vector((nvtProj - fVect).dotProduct(Z),
                                (nvtProj - fVect).dotProduct(Y), 0)
            pvtReduced = vector((pvtProj - fVect).dotProduct(Z),
                                (pvtProj - fVect).dotProduct(Y), 0)
            nfvReduced = vector((projOfN - fVect).dotProduct(Z),
                                (projOfN - fVect).dotProduct(Y), 0)
            pfvReduced = vector((projOfP - fVect).dotProduct(Z),
                                (projOfP - fVect).dotProduct(Y), 0)
            normal = (pvtReduced - nvtReduced).twoDNormal()
            s1 = normal.dotProduct(pfvReduced - pvtReduced)
            if s1 * normal.dotProduct(nfvReduced - pvtReduced) >= 0:
                if s1 > 0:
                    sil.segments.append((pvtReduced.coords2D(),
                                         nvtReduced.coords2D()))
                else:
                    sil.segments.append((nvtReduced.coords2D(),
                                         pvtReduced.coords2D()))
                sil.edgesTaken.append(e)
        return sil

    def visualHull(sils, length):
        """ silhouette list * float -> polyhedron
        Returns the visual hull corresponding to the visual hull
            calculated with the silhouettes <sils>."""
        result = sils.pop(0).cone(length)
        assert result.pnFacesInPoly()
        i = 0
        for s in sils:
            # print(i)
            assert result.pnFacesInPoly()
            result = result.intersection(s.cone(length), True)
            # result.plot()
            i += 1
        return result

    def facesInVertices(self):
        return all(all(v in self.vertices for v in f.vertices) for f in self.faces)

    def pnFacesInPoly(self):
        return all(e.pFace in [None] + self.faces and e.nFace in [None] + self.faces for e in self.edges)

    def nonDoubleVertices(self):
        return (all(all(v1.dist(v2) > COMPARISON_EPSILON or v1 == v2 for v2 in self.vertices) for v1 in self.vertices) and
                all(self.vertices.count(v) == 1 for v in self.vertices))


class intersector:
    """class representing an intersection between an edge and a face"""

    def __init__(self, v, ne, pe, f, faces=[],
                 sTraced=False, iTraced=False,
                 nV=None, nF=None):
        """ vertex * edge * edge * face * face list * bool * bool * vertex * face ->
            intersector
        <pe> and <nf> are the two halves of the edge
            intersecting with the face <f>,
            <v> is the vertex at the point of intersection"""
        v.intersector = self
        assert v.intersector is not None
        self.v = v
        self.pe = pe
        self.ne = ne
        self.f = f
        self.sTraced = sTraced
        self.iTraced = iTraced
        self.nV = nV
        self.nF = nF
        self.faces = faces
        self.adjacents = (None, None)

    def __str__(self):
        return ("""Intersection between edges {} and {} and face {}.
                  \n This intersection is {}sTraced and is {}iTraced"""
                .format(self.pe, self.ne, self.f,
                        "" if self.sTraced else "not ",
                        "" if self.iTraced else "not "))

    def isSTraced(self):
        """intersector -> bool
        Returns wether the intersector has been sTraced."""
        try:
            return self.sTraced
        except AttributeError:
            self.sTraced = False
            return False

    def isITraced(self):
        """intersector -> bool
        Returns wether the intersector has been iTraced."""
        try:
            return self.iTraced
        except AttributeError:
            self.iTraced = False
            return False


class polyIntersection:
    """class used for polyhedron intersection"""

    def __init__(self, poly1, poly2, inter1, inter2, result=None):
        """ polyhedron * polyhedron * intersector list * intersector list *
                polyhedron -> polyIntersection
        <poly1> and <poly2> are the two polyhedron intersecting,
            <intersectors> is the list of their intersectors.
        All elements from <inter1> should have vertices in <poly1>
            and those from <inter2> in <poly2>."""
        self.poly1 = poly1
        self.poly2 = poly2
        self.inter1 = inter1
        self.inter2 = inter2
        if result is not None:
            self.result = result
        else:
            self.result = polyhedron([], [], [])

    def getIntersector(self, v):
        """ polyIntersection * vertex -> intersector/None
        Returns (if it exists) the intersector corresponding to <v>.
        If there is none to be found, returns False"""
        try:
            return v.intersector
        except AttributeError:
            return None

    def getIntersectorList(self, l):
        """ polyIntersection * vertex list -> intersector list
        Returns the list of intersectors corresponding to elements from <l>."""
        return [self.getIntersector(v) for v in l]

    def isIntersection(self, v):
        """ polyIntersection * vertex -> bool
        Returns wether v is a point of intersection."""
        return (any(inter.v == v for inter in self.inter1) or
                any(inter.v == v for inter in self.inter2))

    def addIntersection(self, el1, el2):
        """ face/edge * face/edge -> edge
        Alters <self> to account for the intersection between <el1> and <el2>.
        The intersector object is added to <self>.intersectors,
            and the piercing edge is split.
        <el1> and <el2> must belong to <self>.poly1 and <self>.poly2,
            respectively, and there must be an edge and a face."""
        if type(el1) == edge:
            if type(el2) == face:
                c = el1.faceIntersection(el2)
                if c is not False:
                    assert all(all(v in self.poly1.vertices for v in f.vertices) for f in self.poly1.faces)
                    newV, newEdge = self.poly1.edgeSplit(el1, c)
                    assert all(all(v in self.poly1.vertices for v in f.vertices) for f in self.poly1.faces)
                    self.inter1.append(intersector(newV, el1, newEdge, el2,
                                                   [el1.nFace,
                                                    el1.pFace,
                                                    el2]))
                    return newEdge
                    # # print('intersected')
                el1.markIntersectedWith(el2)
            else:
                raise TypeError('The arguments given must '
                                'be an edge and a face')
        elif type(el1) == face:
            if type(el2) == edge:
                c = el2.faceIntersection(el1)
                if c is not False:
                    newV, newEdge = (self.poly2.edgeSplit(el2, c))
                    self.inter2.append(intersector(newV, el2, newEdge, el1,
                                                   [el2.nFace,
                                                    el2.pFace,
                                                    el1]))
                    return newEdge
                    # # print('intersected')
                el2.markIntersectedWith(el1)
            else:
                raise TypeError('The arguments given must '
                                'be an edge and a face')
        else:
            raise TypeError('The arguments given must be an edge and a face')

    def nextIntersectors(self, inter):
        """ polyIntersection * intersector -> intersector * intersector
        Returns the next intersectors in a surface intersection loop.
        Corresponds approximately to Baumgart's NEXTPV routine."""
        assert self.poly1.pnFacesInPoly() and self.poly2.pnFacesInPoly()
        otherFInters = self.getIntersectorList(inter.f.vertices)
        # First intersector
        pInters = self.getIntersectorList(inter.pe.pFace.vertices)
        otherI1 = next(filter(lambda x: x is not None and x.f == inter.f and x != inter,
                              pInters),
                       None)
        if otherI1 is None:
            # The pFace does not intersect inter.f a second time,
            # looking for a place where inter.f intersects the pFace
            otherI1 = next(filter(lambda x: x is not None and x.f == inter.pe.pFace,
                                  otherFInters),
                           None)
            if otherI1 is None:
                # polyhedron(inter.f.vertices + inter.pe.pFace.vertices,
                #            inter.pe.pFace.edges() + inter.f.edges(),
                #            [inter.f, inter.pe.pFace]).plot(True, col=('none', 'k', 'r'))
                # # print(inter.f, '\n\n', inter.pe.pFace)
                assert all(v in self.poly1.vertices for v in inter.f.vertices) or all(v in self.poly2.vertices for v in inter.f.vertices)
                assert self.poly1.facesInVertices() and self.poly2.facesInVertices()
                assert self.poly1.pnFacesInPoly() and self.poly2.pnFacesInPoly()
                assert self.poly1.nonDoubleVertices() and self.poly2.nonDoubleVertices()
                # # print('\n', [(i, i.f) for i in filter(lambda x: x is not False, pInters)])
                # # print([sum(min(v.dist(v2) for v2 in inter.f.vertices) for v in i.f.vertices) for i in filter(lambda x: x is not False, pInters)])
                # # print('\n', [(i, i.f) for i in filter(lambda x: x is not False, otherFInters)])
                # # print([sum(min(v.dist(v2) for v2 in inter.pe.pFace.vertices) for v in i.f.vertices) for i in filter(lambda x: x is not False, otherFInters)])
                raise ValueError('No intersector found')
        # Second intersector
        nInters = self.getIntersectorList(inter.pe.nFace.vertices)
        otherI2 = next(filter(lambda x: x is not None and x.f == inter.f and x != inter,
                              nInters),
                       None)
        if otherI2 is None:
            # The nFace does not intersect inter.f a second time,
            # looking for a place where inter.f intersects the nFace
            otherI2 = next(filter(lambda x: x is not None and x.f == inter.pe.nFace,
                                  otherFInters),
                           None)
            if otherI2 is None:
                polyhedron(inter.f.vertices + inter.pe.nFace.vertices,
                           inter.pe.nFace.edges() + inter.f.edges(),
                           [inter.f, inter.pe.nFace]).plot(True, col=('none', 'k', 'r'))
                # # print(inter.f, inter.pe.pFace)
                raise ValueError('No intersector found')
        inter.adjacents = (otherI1, otherI2)
        return (otherI1, otherI2)

    def fromPolyhedron(poly1, poly2):
        """ polyhedron * polyhedron -> polyIntersection
        Builds the polyIntersection object for poly1 and poly2."""
        polyInt = polyIntersection(poly1, poly2, [], [])
        originalEdges1 = copy.copy(poly1.edges)
        originalEdges2 = copy.copy(poly2.edges)
        assert all(all(v in poly1.vertices for v in f.vertices) for f in poly1.faces)
        assert all(f.hasNoDoubleVertices() for f in polyInt.poly1.faces)
        assert all(f.hasNoDoubleVertices() for f in polyInt.poly2.faces)
        assert all(all(v in polyInt.poly1.vertices for v in f.vertices) for f in polyInt.poly1.faces)
        assert all(all(v in polyInt.poly2.vertices for v in f.vertices) for f in polyInt.poly2.faces)
        while len(originalEdges1) + len(originalEdges2) > 0:
            newEdges1 = []
            newEdges2 = []
            for e1 in originalEdges1:  # WARNING: The list of edges can change, things can get messy.
                for f2 in poly2.faces:
                    if not e1.hasIntersectedWith(f2):
                        assert all(all(v in poly1.vertices for v in f.vertices) for f in poly1.faces)
                        assert all(f.hasNoDoubleVertices() for f in polyInt.poly1.faces)
                        assert all(f.hasNoDoubleVertices() for f in polyInt.poly2.faces)
                        try:
                            assert polyInt.nonDoubleIntersectors()
                        except:
                            pass
                            # # print('second test:', next((str(x[0]), str(x[1])) for x in itertools.product(polyInt.inter1 + polyInt.inter2, polyInt.inter1 + polyInt.inter2) if x[0] != x[1] and x[0].v == x[1].v and x[0].f == x[1].f))
                            # print(next((str(x[0]), str(x[1])) for x in itertools.product(polyInt.inter1 + polyInt.inter2, polyInt.inter1 + polyInt.inter2) if x[0] != x[1] and x[0].v == x[1].v))
                            # assert False
                        ne = polyInt.addIntersection(e1, f2)
                        if ne is not None:
                            ne.markIntersectedWith(f2)
                            newEdges1.append(ne)
            for e2 in originalEdges2:
                for f1 in poly1.faces:
                    if not e2.hasIntersectedWith(f1):
                        ne = polyInt.addIntersection(f1, e2)
                        if ne is not None:
                            ne.markIntersectedWith(f1)
                            newEdges1.append(ne)
            originalEdges1 = newEdges1
            originalEdges2 = newEdges2
            assert all(all(v in polyInt.poly1.vertices for v in f.vertices) for f in polyInt.poly1.faces)
            assert all(all(v in polyInt.poly2.vertices for v in f.vertices) for f in polyInt.poly2.faces)
        return polyInt

    def buildEdges(self):
        """ polyIntersection -> None
        Builds the edges of the result."""
        assert self.poly1.pnFacesInPoly() and self.poly2.pnFacesInPoly()
        assert self.nonDoubleIntersectors()
        # print('2.1')
        for i in self.inter1:
            assert self.nonDoubleIntersectors()
            # print('test')
            assert all(all(v in self.poly1.vertices for v in f.vertices) for f in self.poly1.faces)
            if not i.isSTraced():
                self.traceSurface(i)
        # print('2.2')
        for i in self.inter2:
            if not i.isSTraced():
                self.traceSurface(i)
        for i in self.inter1:
            if not i.isITraced():
                self.traceInterior(i)
        for i in self.inter2:
            if not i.isITraced():
                self.traceInterior(i)

    def traceSurface(self, inter):
        """ polyIntersection * intersector -> ()
        Traces the surface loop that <inter> belongs to,
        and adds it to <self>.result.
        Marks the intersector that have been processed."""
        assert self.nonDoubleIntersectors()
        assert self.poly1.pnFacesInPoly() and self.poly2.pnFacesInPoly()
        assert all(f.hasNoDoubleVertices() for f in self.poly1.faces)
        assert all(f.hasNoDoubleVertices() for f in self.poly2.faces)
        assert self.nonDoubleVertices()
        i = self.nextIntersectors(inter)[0]
        assert i.v != inter.v
        assert i.v.dist(inter.v) != 0
        lastI = inter
        self.result.addEdge(self.newVertex(lastI), self.newVertex(i))
        i.sTraced = True
        while i != inter:
            nextInts = self.nextIntersectors(i)
            if nextInts[0] != lastI:
                lastI = i
                i = nextInts[0]
            else:
                lastI = i
                i = nextInts[1]
            self.result.addEdge(self.newVertex(lastI), self.newVertex(i))
            i.sTraced = True

    def newVertex(self, inter):
        """ polyIntersection * intersector/vertex -> vertex
        Returns the new vertex for <inter> (and creates it if needed)."""
        try:
            if inter.nV is None:
                raise AttributeError
            return inter.nV
        except AttributeError:
            vert = inter.v if type(inter) is intersector else inter
            inter.nV = self.result.addVertex(vert.x, vert.y, vert.z)
            assert inter.nV is not None
            return inter.nV

    def createIntersectorVertices(self):
        """ polyIntersection -> None
        Adds the vertices corresponding to the intersectors to the result."""
        for i in self.inter1 + self.inter2:
            self.result.addVertex(self.newVertex(i))

    def traceInterior(self, inter):
        """ polyIntersection * intersector -> ()
        Traces the interior graph that <inter> belongs to,
            and adds it to <self>.result.
        Marks the intersector that have been processed."""
        v1 = vector(inter.pe.other(inter.v)) - vector(inter.v)
        normal = inter.f.normalVect()
        if v1.dotProduct(normal) < 0:
            vert = inter.pe.other(inter.v)
        else:
            vert = inter.ne.other(inter.v)
        nI = self.getIntersector(vert)
        if nI is not None:
            self.result.addEdge(self.newVertex(inter), self.newVertex(nI))
            nI.iTraced = True
        else:
            if inter.v in self.poly1.vertices:
                poly = self.poly1
            else:
                poly = self.poly2
            newVert = self.newVertex(vert)
            vert.markInt()
            newVert.markInt()
            queue = [(vert, newVert)]
            done = []
            while len(queue) > 0:
                vert, newVert = queue.pop()
                self.result.addVertex(newVert)
                for e in poly.edges:
                    if vert in [e.pvt, e.nvt]:
                        otherV = e.other(vert)
                        otherI = self.getIntersector(otherV)
                        if otherI is not None:
                            otherNV = self.newVertex(otherI)
                            self.result.addEdge(newVert, otherNV)
                            otherI.iTraced = True
                        elif otherV not in done:
                            otherNV = self.newVertex(otherV)
                            otherNV.markInt()
                            otherV.markInt()
                            self.result.addEdge(newVert, otherNV)
                            queue.append((otherV, otherNV))
                done.append(vert)
        inter.iTraced = True

    def faceTrace(self, pi, fi, vi):
        """ int * int * int -> None
        Builds the face of the result that self.poly<pi>.faces[fi].vertices[vi]
            belongs to.
        This vertex must be an intersection."""
        if pi == 1:
            face = self.poly1.faces[fi]
        elif pi == 2:
            face = self.poly2.faces[fi]
        else:
            raise ValueError('pi must be 1 or 2')
        n = len(face.vertices)
        initialV = face.vertices[vi]
        newFace = [initialV]
        initialI = self.getIntersector(initialV)
        assert initialI is not None  # This checks that the vertex was actually an intersection
        # This part determines which side of the face is the way to go
        if self.getIntersector(face.vertices[(vi + 1) % n]) is not None or face.vertices[(vi + 1) % n].isInterior():
            step = 1
        elif self.getIntersector(face.vertices[(vi - 1) % n]) is not None or face.vertices[(vi - 1) % n].isInterior():
            step = -1
        else:
            raise ValueError('The vertex given has no adjacent interior or intersector vertex')
        i = (vi + step) % n
        force = True
        while force or i != (vi + step) % n:
            force = False
            # Tracing the interior
            while (face.vertices[i]).isInterior():
                newFace.append(self.newVertex(face.vertices[i]))
                face.vertices[i].markFTraced(face)
                # print('fTracing', i)
                i += step
            # Tracing the surface
            newFace.append(self.newVertex(face.vertices[i]))
            face.vertices[i].markFTraced(face)
            lastI = self.getIntersector(face.vertices[i])
            potentialInters = self.nextIntersectors(lastI)
            if abs((vector(potentialInters[0].v) - vector(lastI.v))
               .dotProduct(face.normalVect())) < COMPARISON_EPSILON:
               inter = potentialInters[0]
            elif abs((vector(potentialInters[1].v) - vector(lastI.v))
               .dotProduct(face.normalVect())) < COMPARISON_EPSILON:
               inter = potentialInters[1]
            else:
                raise Exception('The nextIntersectors are weird (1)')
            while inter.v not in face.vertices:
                nextInters = self.nextIntersectors(inter)
                if nextInters[0] != lastI:
                    lastI = inter
                    inter = nextInters[0]
                elif nextInters[1] != lastI:
                    lastI = inter
                    inter = nextInters[1]
                else:
                    raise Exception('The nextIntersectors are weird (2)')
                newFace.append(lastI.v)
            inter.v.markFTraced(face)
            i = (face.vertices.index(inter.v) + step) % n
        if step == -1:
            newFace.reverse()
        self.result.addFace(newFace, True)

    def buildFaces(self):
        """ polyIntersection -> None
        Builds the faces of the result of the intersection."""
        for fi in range(len(self.poly1.faces)):
            f = self.poly1.faces[fi]
            for vi in range(len(f.vertices)):
                v = f.vertices[vi]
                if self.getIntersector(v) is not None and not v.isMarkedFor(f):
                    self.faceTrace(1, fi, vi)
            v = f.vertices[0]
            if v.isInterior() and not v.isMarkedFor(f):
                interiorFaceTrace(f)
        for fi in range(len(self.poly2.faces)):
            f = self.poly2.faces[fi]
            for vi in range(len(f.vertices)):
                v = f.vertices[vi]
                if self.getIntersector(v) is not None and not v.isMarkedFor(f) :
                    self.faceTrace(2, fi, vi)
            v = f.vertices[0]
            if v.isInterior() and not v.isMarkedFor(f):
                interiorFaceTrace(f)

    def interiorFaceTrace(self, f):
        """ polyIntersection * face -> None
        Adds the new face corresponding to <f> in the result.
        <f> must be completely interior."""
        assert all(v.isInterior() for v in f.vertices)
        self.result.addFace([self.newVertex(v) for v in f.vertices])

    def nonDoubleVertices(self):
        return self.poly1.nonDoubleVertices() and self.poly2.nonDoubleVertices()

    def nonDoubleIntersectors(self):
        return True
        # return all(all(i1 == i2 or i1.v != i2.v for i1 in self.inter1 + self.inter2) for i2 in self.inter1 + self.inter2)


class vector:
    """Simple three-dimensional vectors"""

    def __init__(self, x=0, y=0, z=0):
        if type(x) is tuple:
            self.x = x[0]
            self.y = x[1]
            self.z = x[2]
        elif type(x) is vertex:
            self.x = x.x
            self.y = x.y
            self.z = x.z
        elif type(x) is float or type(x) is int:
            self.x = x
            self.y = y
            self.z = z
        elif type(x) is edge:
            self.x = x.pvt.x - x.nvt.x
            self.y = x.pvt.y - x.nvt.y
            self.z = x.pvt.z - x.nvt.z
        else:
            raise ValueError('Invalid type given to the vector constructor')

    def __str__(self):
        return "({}, {}, {})".format(self.x, self.y, self.z)

    def __mul__(self, el2):
        return el2.__rmul__(self)

    def __rmul__(self, el2):
        """ float / vector -> vector
        If <el2> is scalar,
            returns <v> with each coordinate multiplied by <el2>.
        If <el2> is a vector,
            returns the cross product of the vectors <el2> and <self>."""
        if type(el2) is float or type(el2) is int:
            return vector(el2 * self.x, el2 * self.y, el2 * self.z)
        elif type(el2) is vector:
            return vector(el2.y * self.z - el2.z * self.y,
                          el2.z * self.x - el2.x * self.z,
                          el2.x * self.y - el2.y * self.x)
        else:
            raise TypeError('Cannot multiply a vector with something'
                            'that is neither a vector, a float or an int')

    def __add__(self, v):
        """ vector * vector -> vector
        Returns the sum of the two vectors"""
        return vector(self.x + v.x, self.y + v.y, self.z + v.z)

    def __sub__(self, v):
        """ vector * vector -> vector
        Returns the difference of the two vectors"""
        return self + (-1) * v

    def norm(self):
        """ vector -> float
        Returns the norm of <self>."""
        return math.sqrt(self.dotProduct(self))

    def coords(self):
        """ vector -> float**3
        Returns the tuple of coordinates of <self>."""
        return (self.x, self.y, self.z)

    def coords2D(self):
        """ vector -> float**2
        Returns the tuple of the first two coordinates of <self>."""
        return (self.x, self.y)

    def dotProduct(self, v):
        """ vector * vector -> float
        Returns the dot product of the two vectors."""
        return self.x * v.x + self.y * v.y + self.z * v.z

    def normalize(self):
        """ vector -> vector"""
        try:
            return (1 / self.norm()) * self
        except ZeroDivisionError:
            raise ValueError('Cannot normalize a null vector')

    def xCrossProd(self, other):
        """ vector * vector -> float
        Returns the x component of the cross product of self and other."""
        return other.y * self.z - other.z * self.y

    def twoDNormal(self):
        """ vector -> vector
        Returns a vector normal to <self>."""
        return vector((-1) * self.y, self.x, 0)


class silhouette:
    """Class representing the silhouette obtained from a photography."""

    def __init__(self, focalPoint, focalDist, angles, segments):
        """ (float * float * float) * float * (float  * float * float) *
            (float * float)**2 list -> silhouette"""
        self.focalPoint = focalPoint
        self.focalDist = focalDist
        self.angles = angles
        self.segments = segments

    def cone(self, length):
        """ silhouette * float -> polyhedron
        Returns the cone polyhedron corresponding to the silhouette."""
        X, Y, Z = baseFromAngles(self.angles)
        result = polyhedron([], [], [])
        assert result.pnFacesInPoly()
        f = result.addVertex(self.focalPoint)
        for s in self.segments:
            vect0 = (self.focalDist * X +
                     s[0][0] * Z +
                     s[0][1] * Y).normalize()
            vect1 = (self.focalDist * X +
                     s[1][0] * Z +
                     s[1][1] * Y).normalize()
            assert result.pnFacesInPoly()
            vert0 = result.addVertex(vector(self.focalPoint) +
                                     length * vect0)
            assert result.pnFacesInPoly()
            vert1 = result.addVertex(vector(self.focalPoint) +
                                     length * vect1)
            assert result.pnFacesInPoly()
            try:
                result.addFace([f, vert1, vert0])
            except ValueError:
                result.plot()
                p = polyhedron([], [], [])
                p.addFace([f, vert1, vert0])
                p.union(result).plot()
                self.plot()
                assert False
            assert result.pnFacesInPoly()
        # WARNING : the cone is not closed at its top.
        # I'm not sure if this can cause issues.
        return result

    def plot(self):
        """ silhouette -> None
        Plots the silhouette using matplotlib"""
        done = []
        patches = []

        def chain(p):
            done.append(p)
            yield p
            last = p
            finished = False
            while not finished:
                seg0 = next(filter(lambda x: (x[0] not in done and
                                              x[1] == last),
                                   self.segments),
                            None)
                seg1 = next(filter(lambda x: (x[1] not in done and
                                              x[0] == last),
                                   self.segments),
                            None)
                if seg0 is not None:
                    done.append(seg0[0])
                    yield seg0[0]
                    last = seg0[0]
                elif seg1 is not None:
                    done.append(seg1[1])
                    yield seg1[1]
                    last = seg1[1]
                else:
                    finished = True
        p0 = next(filter(lambda x: x not in done,
                         (s[0] for s in self.segments)),
                  next(filter(lambda x: x not in done,
                              (s[1] for s in self.segments)),
                       None))
        assert type(p0) is tuple
        while p0 is not None:
            gen = chain(p0)
            a = np.array([[p[0], p[1]] for p in gen])
            patches.append(mplPoly(a, closed=True))
            p0 = next(filter(lambda x: x not in done,
                             (s[0] for s in self.segments)),
                      next(filter(lambda x: x not in done,
                                  (s[1] for s in self.segments)),
                           None))
        fig, ax = plt.subplots()
        ax.set_xlim([min(min(p.get_xy()[:, 0]) for p in patches),
                     max(max(p.get_xy()[:, 0]) for p in patches)])
        ax.set_ylim([min(min(p.get_xy()[:, 1]) for p in patches),
                     max(max(p.get_xy()[:, 1]) for p in patches)])
        pCol = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
        colors = 100 * np.random.rand(len(patches))
        pCol.set_array(np.array(colors))
        ax.add_collection(pCol)
        plt.show()

    def __str__(self):
        return str(self.segments)


def splitList(l, i1, i2):
    """ a' list * a' * a' -> a' list * a' list
    Splits <l> into two lists, the first going from <i1> to <i2>,
        the second going from <i2> to <i1> (both loop through the end of <l>).
    The lists returned both have <l>[<i1>] and <l>[<i2>]."""
    n = len(l)
    L1 = []
    L2 = []
    j = 0
    # Creating list 1
    while (i1 + j) % n != i2:
        L1.append(l[(i1 + j) % n])
        j += 1
    L1.append(l[i2])
    # Creating list 2
    j = 0
    while (i2 + j) % n != i1:
        L2.append(l[(i2 + j) % n])
        j += 1
    L2.append(l[i1])
    return (L1, L2)


def planeLineIntersect(p1, p2, equ):
    """ float**3 * float**3 * float**4 -> float**3
    Returns the intersection of the line that passes through the points
        of coordinates p1 and p2 and the plane of equation
        <equ>[0]*x + <equ>[1]*y + <equ>[2]*z = <equ>[3]"""
    n = vector(equ[0], equ[1], equ[2])
    v1, v2 = vector(p1), vector(p2)
    t = (equ[3] - n.dotProduct(v2)) / (n.dotProduct(v1 - v2))
    return (t * v1 + (1 - t) * v2).coords()


def onBase():
    """ None -> vector**3
    Returns an orthonormal base of R**3."""
    return (vector(1, 0, 0), vector(0, 1, 0), vector(0, 0, 1)) 


def baseFromAngles(angles, originalBase=onBase()):
    """ float**3 -> vector**3
    Returns the base that has the angles of rotation <angles>
        relative to originalBase.
    The segments are defined in the coordinates system that has :
       self.focalPoint + focalDist * X as origin
       Z and Y as basis;"""
    x, y, z = originalBase
    yaw, roll, pitch = angles
    P = math.cos(yaw) * x + math.sin(yaw) * y
    x2 = math.cos(pitch) * P + math.sin(pitch) * z
    y2 = (z * x2).normalize()
    z2 = x2 * y2
    y3 = math.cos(roll) * y2 - math.sin(roll) * z2
    # N = math.cos(yaw) * x + math.sin(yaw) * y
    # A = (N * z).normalize()
    # Z = math.cos(roll) * z + math.sin(roll) * A
    # B = (Z * N).normalize()
    # X = math.cos(pitch) * N + math.sin(pitch) * B
    # plotVects([x, z, N, A, B, X])
    # Y = Z * X
    return (x2, y3, x2 * y3)


def simpleSphere(precision):
    """ int -> polyhedron
    Returns an approximation of sphere with 2*<precision>**2 faces."""
    b = polyhedron([vertex(0, 0, 1)], [], [])
    rot1 = vector(math.pi / precision, 0, 0)
    for i in range(2 * precision + 1):
        b.rotateSweep(rot1, (0, 0, 0))
    rot2 = vector(0, 0, math.pi / precision)
    for j in range(precision):
        b.rotateSweep(rot2, (0, 0, 0))
    return b


def simpleTetrahedron():
    v1 = vertex(0, 0, 0)
    v2 = vertex(1, 0, 0)
    v3 = vertex(0, 1, 0)
    v4 = vertex(0, 0, 1)
    tetra = polyhedron([v1, v2, v3, v4], [], [])
    tetra.addEdge(v1, v2)
    tetra.addEdge(v1, v3)
    tetra.addEdge(v1, v4)
    tetra.addEdge(v2, v3)
    tetra.addEdge(v3, v4)
    tetra.addEdge(v4, v2)
    tetra.addFace([v1, v2, v3])
    tetra.addFace([v1, v3, v4])
    tetra.addFace([v1, v2, v4])
    tetra.addFace([v2, v3, v4])
    return tetra


def simpleCube():
    v1 = vertex(0, 0, 0)
    v2 = vertex(1, 0, 0)
    v3 = vertex(0, 1, 0)
    v4 = vertex(0, 0, 1)
    v5 = vertex(1, 1, 0)
    v6 = vertex(1, 0, 1)
    v7 = vertex(0, 1, 1)
    v8 = vertex(1, 1, 1)
    cube = polyhedron([v1, v2, v3, v4, v5, v6, v7, v8], [], [])
    cube.addEdge(v1, v2)
    cube.addEdge(v1, v3)
    cube.addEdge(v1, v4)
    cube.addEdge(v2, v5)
    cube.addEdge(v3, v5)
    cube.addEdge(v4, v6)
    cube.addEdge(v2, v6)
    cube.addEdge(v4, v6)
    cube.addEdge(v3, v7)
    cube.addEdge(v4, v7)
    cube.addEdge(v5, v8)
    cube.addEdge(v6, v8)
    cube.addEdge(v7, v8)
    cube.addFace([v1, v3, v5, v2])
    cube.addFace([v1, v2, v6, v4])
    cube.addFace([v1, v4, v7, v3])
    cube.addFace([v2, v5, v8, v6])
    cube.addFace([v3, v7, v8, v5])
    cube.addFace([v4, v6, v8, v7])
    return cube


def rotateFunction(rotCenter, rotVect):
    """ (float * float * float) * vector
        -> ((float * float * float) -> (float * float * float))
        Returns the function for rotating a set of coordinates."""
    rotn = (1 / rotVect.norm()) * rotVect

    def rot(coords):
        oldVect = vector(coords) - vector(rotCenter)
        newVect = (math.cos(rotVect.norm()) * oldVect +
                   math.sin(rotVect.norm()) * rotn * oldVect +
                   ((1 - math.cos(rotVect.norm())) *
                    (rotn.dotProduct(oldVect)) * rotn))
        return (rotCenter[0] + newVect.x,
                rotCenter[1] + newVect.y,
                rotCenter[2] + newVect.z)
    return rot


def simplex(v1, v2, v3):
    """ vector * vector * vector -> polyhedron
    Returns a simplex made with the three vector."""
    res = polyhedron([], [], [])
    p0 = res.addVertex(vertex(0, 0, 0))
    p1 = res.addVertex(vertex(v1))
    p2 = res.addVertex(vertex(v2))
    p3 = res.addVertex(vertex(v3))
    p12 = res.addVertex(vertex(v1 + v2))
    p13 = res.addVertex(vertex(v1 + v3))
    p23 = res.addVertex(vertex(v2 + v3))
    p123 = res.addVertex(vertex(v1 + v2 + v3))
    res.addEdge(p0, p1)
    res.addEdge(p0, p2)
    res.addEdge(p0, p3)
    res.addEdge(p1, p12)
    res.addEdge(p1, p13)
    res.addEdge(p2, p12)
    res.addEdge(p2, p23)
    res.addEdge(p3, p13)
    res.addEdge(p3, p23)
    res.addEdge(p1, p123)
    res.addEdge(p2, p123)
    res.addEdge(p3, p123)
    res.addFace([p0, p2, p12, p1])
    res.addFace([p0, p1, p13, p3])
    res.addFace([p0, p3, p23, p2])
    res.addFace([p1, p12, p123, p13])
    res.addFace([p2, p23, p123, p12])
    res.addFace([p3, p13, p123, p23])
    return res


def plotVects(vList, colors='k'):
    """ vector list -> None
    Plots a list of vectors."""
    polyhedron([vertex(v) for v in vList] + [vertex(0, 0, 0)], [edge(v) for v in vList], []).plot(plotEdges=True)

# WARNING: not working ?
def insert(l, i, el):
    l.append(l[-1])
    for j in range(len(l)-1, i+1, -1):
        l[j] = l[j-1]
    l[i] = el
    return l
