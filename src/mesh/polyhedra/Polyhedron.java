package mesh.polyhedra;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Vector3d;

import mesh.Edge;
import mesh.Face;
import mesh.Mesh;
import mesh.struct.OrderedVertexToAdjacentEdge;
import mesh.struct.OrderedVertexToAdjacentFace;

/**
 * A class for generic convex closed polyhedra meshes. The class contains
 * several operations for creating new polyhedra from existing ones. The
 * operation names and behaviors are based on Conway's polyhedron notation.
 * Each such operation creates an entirely new mesh (the original remains
 * unchanged by the operations).
 * 
 * @author Brian Yao
 */
public class Polyhedron extends Mesh {
	
	/**
	 * Computes the dual polyhedron of this polyhedron. The dual polyhedron
	 * has the "opposite" topology, in that it has one vertex for each face of
	 * the original polyhedron, and one face for each vertex of the original
	 * polyhedron.
	 * 
	 * @return The dual polyhedron.
	 */
	public Polyhedron dual() {
		Polyhedron dualPolyhedron = new Polyhedron();
		
		// Create new vertices, one at the centroid of each face
		Map<Face, Integer> newVertices = new HashMap<>();
		int vertexIndex = 0;
		for (Face face : this.getFaces()) {
			Vector3d cent = face.centroid();
			cent.scale(1.0 / cent.lengthSquared());
			dualPolyhedron.addVertexPosition(cent);
			newVertices.put(face, vertexIndex++);
		}
		
		// Construct new faces
		OrderedVertexToAdjacentFace ovtaf = new OrderedVertexToAdjacentFace(this);
		for (int i = 0 ; i < numVertexPositions() ; i++) {
			if (ovtaf.getAdjacentFaces(i) != null) {
				List<Face> adjacentFaces = ovtaf.getAdjacentFaces(i);
				Face newFace = new Face(adjacentFaces.size());
				for (int j = 0 ; j < adjacentFaces.size() ; j++) {
					Face jthFace = adjacentFaces.get(j);
					int newVertexIndex = newVertices.get(jthFace);
					newFace.setVertexPosition(j, newVertexIndex);
				}
				dualPolyhedron.addFace(newFace);
			}
		}
		
		dualPolyhedron.setVertexNormalsToFaceNormals();
		return dualPolyhedron;
	}
	
	/**
	 * Computes the rectification of this polyhedron. This involves creating
	 * new vertices at the midpoints of the edges of the original polyhedron,
	 * and discarding the vertices of the original (they become faces). "Ambo"
	 * refers to the name of the operation under Conway's notation.
	 * 
	 * @return The ambo polyhedron.
	 */
	public Polyhedron ambo() {
		Polyhedron amboPolyhedron = new Polyhedron();
		
		// Create new vertices, one at the midpoint of each edge
		Map<Edge, Integer> newVertices = new HashMap<>();
		int vertexIndex = 0;
		for (Face face : this.getFaces()) {
			// Create a new face for each face on this polyhedron
			Face newFace = new Face(face.numVertices());
			Edge[] edges = face.getEdges();
			for (int i = 0 ; i < face.numVertices() ; i++) {
				Edge edge = edges[i];
				if (newVertices.get(edge) == null) {
					Vector3d edgeMidpt = edge.midpoint();
					amboPolyhedron.addVertexPosition(edgeMidpt);
					newVertices.put(edge, vertexIndex++);
				}
				newFace.setVertexPosition(i, newVertices.get(edge));
			}
			amboPolyhedron.addFace(newFace);
		}
		
		// Construct new faces in place of vertices of original
		OrderedVertexToAdjacentEdge ovtae = new OrderedVertexToAdjacentEdge(this);
		for (int i = 0 ; i < numVertexPositions() ; i++) {
			List<Edge> adjacentEdges = ovtae.getAdjacentEdges(i);
			Face newVertexFace = new Face(adjacentEdges.size());
			for (int j = 0 ; j < adjacentEdges.size() ; j++) {
				newVertexFace.setVertexPosition(j, newVertices.get(adjacentEdges.get(j)));
			}
			amboPolyhedron.addFace(newVertexFace);
		}
		
		amboPolyhedron.setVertexNormalsToFaceNormals();
		return amboPolyhedron;
	}
	
	/**
	 * Computes the "join" polyhedron of this polyhedron. Each face of the
	 * original polyhedron becomes a raised pyramid, such that for two adjacent
	 * pyramids, the faces of the pyramids which share an edge are coplanar. In
	 * this way, we obtain a polyhedron made entirely of quadrilateral faces.
	 * 
	 * @return The join polyhedron.
	 */
	public Polyhedron join() {
		return this.ambo().dual();
	}
	
}
