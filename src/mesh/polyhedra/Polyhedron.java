package mesh.polyhedra;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Vector3d;

import util.Canonicalize;
import mesh.Edge;
import mesh.Face;
import mesh.Mesh;
import mesh.struct.EdgeToAdjacentFace;
import mesh.struct.OrderedVertexToAdjacentEdge;
import mesh.struct.OrderedVertexToAdjacentFace;

/**
 * A class for generic closed polyhedra meshes. The class contains
 * several operations for creating new polyhedra from existing ones. The
 * operation names and behaviors are based on Conway's polyhedron notation.
 * Each such operation creates an entirely new mesh (the original remains
 * unchanged by the operations).
 * 
 * These operations rely heavily on the convention that faces' vertices are
 * all specified in counter-clockwise order.
 * 
 * @author Brian Yao
 */
public class Polyhedron extends Mesh {
	
	/**
	 * Create an empty polyhedron.
	 */
	public Polyhedron() {
		super();
	}
	
	/**
	 * Copy constructor.
	 * 
	 * @param polyhedron The polyhedron to copy.
	 */
	public Polyhedron(Polyhedron polyhedron) {
		super();
		for (Vector3d vertexPos : vertexPositions) {
			addVertexPosition(new Vector3d(vertexPos));
		}
		for (Vector3d normal : vertexNormals) {
			addVertexNormal(new Vector3d(normal));
		}
		for (Face face : faces) {
			addFace(new Face(face));
		}
	}
	
	/**
	 * Clones a polyhedron. Equivalent to the copy constructor. Modifications
	 * made to the clone are separate from those made to the original.
	 * 
	 * @return A clone of this polyhedron.
	 */
	public Polyhedron clone() {
		Polyhedron clone = new Polyhedron();
		for (Vector3d vertexPos : vertexPositions) {
			clone.addVertexPosition(new Vector3d(vertexPos));
		}
		for (Vector3d normal : vertexNormals) {
			clone.addVertexNormal(new Vector3d(normal));
		}
		for (Face face : faces) {
			clone.addFace(new Face(face));
		}
		return clone;
	}
	
	/**
	 * Canonicalizes this polyhedron. See util.Canonicalize for more details.
	 * 
	 * @param iterations The number of iterations to run the canonicalization.
	 * @return The canonicalized version of this polyhedron.
	 */
	public Polyhedron canonicalize(int iterations) {
		Polyhedron canonicalized = this.clone();
		Canonicalize.adjust(canonicalized, iterations);
		return canonicalized;
	}
	
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
		for (Face face : faces) {
			Vector3d cent = face.vertexAverage();
			cent.scale(1.0 / cent.lengthSquared()); // Necessary for faces to be planar
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
		for (Face face : faces) {
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
	
	/**
	 * Computes the "kis" polyhedron of this polyhedron. Each face of the
	 * original polyhedron becomes a raised pyramid; the returned polyhedron
	 * is not necessarily convex.
	 * 
	 * @return The kis polyhedron
	 */
	public Polyhedron kis() {
		return this.dual().needle();
	}
	
	/**
	 * Computes the truncated polyhedron of this polyhedron. Each vertex is
	 * truncated, leaving behind a polygon face instead.
	 * 
	 * @return The truncated polyhedron.
	 */
	public Polyhedron truncate() {
		return this.needle().dual();
	}
	
	/**
	 * Computes the "needle" polyhedron of this polyhedron. The centroid of
	 * each face becomes a vertex, and for each edge, there are two
	 * triangular faces which cross the original edge; the edge shared by
	 * these two triangular faces goes between the centroids of the two
	 * original faces. The returned polyhedron is not necessarily convex.
	 * 
	 * @return The needle polyhedron.
	 */
	public Polyhedron needle() {
		Polyhedron needlePolyhedron = new Polyhedron();
		needlePolyhedron.addVertexPositions(vertexPositions);
		
		// Construct new vertices, one at the centroid of each face
		Map<Face, Integer> newVertices = new HashMap<>();
		int vertexIndex = vertexPositions.size();
		for (Face face : faces) {
			Vector3d cent = face.centroid();
			needlePolyhedron.addVertexPosition(cent);
			newVertices.put(face, vertexIndex++);
		}
		
		// Construct new faces, two per edge
		EdgeToAdjacentFace etaf = new EdgeToAdjacentFace(this);
		for (Edge edge : etaf.getEdges()) {
			int[] edgeEnds = edge.getEnds();
			Face[] adjacentFaces = etaf.getAdjacentFaces(edge);
			
			// Check orientation: make sure faces we make are specified in
			// CCW order
			for (int i = 0 ; i < adjacentFaces[0].numVertices() ; i++) {
				int v = adjacentFaces[0].getVertexPosition(i);
				if (v == edgeEnds[0]) {
					int vnext = (i + 1) % adjacentFaces[0].numVertices();
					if (adjacentFaces[0].getVertexPosition(vnext) != edgeEnds[1]) {
						// Swap order of adjacent faces
						Face temp = adjacentFaces[0];
						adjacentFaces[0] = adjacentFaces[1];
						adjacentFaces[1] = temp;
					}
					break;
				}
			}
			
			Face topFace = new Face(3);
			Face bottomFace = new Face(3);
			int vface0 = newVertices.get(adjacentFaces[0]);
			int vface1 = newVertices.get(adjacentFaces[1]);
			topFace.setAllVertexPositions(vface0, vface1, edgeEnds[1]);
			bottomFace.setAllVertexPositions(vface1, vface0, edgeEnds[0]);
			needlePolyhedron.addFaces(topFace, bottomFace);
		}
		
		needlePolyhedron.setVertexNormalsToFaceNormals();
		return needlePolyhedron;
	}
	
	/**
	 * Computes the "zip" polyhedron of this polyhedron. Also known as a
	 * bitruncation operation, this is equivalent to the truncation of
	 * the dual polyhedron.
	 * 
	 * @return The zip polyhedron.
	 */
	public Polyhedron zip() {
		return this.kis().dual();
	}
	
	/**
	 * Computes the expanded polyhedron of this polyhedron. Each vertex becomes
	 * a face, and each edge becomes a quadrilateral face. The appearance of
	 * the returned polyhedron looks as if the original polyhedron's faces had
	 * been "pulled apart", and additional faces filled the gaps.
	 * 
	 * @return The expanded polyhedron.
	 */
	public Polyhedron expand() {
		return this.ambo().ambo();
	}
	
	/**
	 * Computes the "ortho" polyhedron of this polyhedron. Each face is divided
	 * into n quadrilateral faces, where n is the number of vertices in the
	 * face.
	 * 
	 * @return The ortho polyhedron.
	 */
	public Polyhedron ortho() {
		return this.expand().dual();
	}
	
	/**
	 * Computes the "gyro" polyhedron of this polyhedron. Each face is divided
	 * into n pentagonal faces, where n is the number of vertices in the face.
	 * 
	 * @return The gyro polyhedron.
	 */
	public Polyhedron gyro() {
		Polyhedron gyroPolyhedron = new Polyhedron();
		for (Vector3d vertexPos : vertexPositions) {
			gyroPolyhedron.addVertexPosition(new Vector3d(vertexPos));
		}
		
		Map<Integer, Map<Integer, int[]>> newEdges = new HashMap<>();
		
		// Generate one vertex per face
		Map<Face, Integer> centroidIndices = new HashMap<>();
		for (Face face : faces) {
			// The centroid of each face becomes a vertex
			Vector3d centroid = face.centroid();
			int centroidIndex = gyroPolyhedron.vertexPositions.size();
			gyroPolyhedron.addVertexPosition(centroid);
			centroidIndices.put(face, centroidIndex);
		}
		
		for (Edge edge : this.getEdges()) {
			// Generate two new vertices per edge
			Vector3d[] endPositions = edge.getEndLocations();
			Vector3d diff = new Vector3d();
			diff.sub(endPositions[1], endPositions[0]);
			diff.scale(1.0 / 3.0);
			
			Vector3d firstNewVertex = new Vector3d();
			firstNewVertex.add(endPositions[0], diff);
			
			Vector3d secondNewVertex = new Vector3d();
			secondNewVertex.add(firstNewVertex, diff);
			
			int firstIndex = gyroPolyhedron.vertexPositions.size();
			int secondIndex = firstIndex + 1;
			gyroPolyhedron.addVertexPositions(firstNewVertex, secondNewVertex);
			
			// Map the existing edge to the new vertices along it
			int[] ends = edge.getEnds();
			if (newEdges.get(ends[0]) == null) {
				newEdges.put(ends[0], new HashMap<Integer, int[]>());
			}
			newEdges.get(ends[0]).put(ends[1], new int[] {firstIndex, secondIndex});
			
			if (newEdges.get(ends[1]) == null) {
				newEdges.put(ends[1], new HashMap<Integer, int[]>());
			}
			newEdges.get(ends[1]).put(ends[0], new int[] {secondIndex, firstIndex});
		}
		
		// Construct pentagonal faces
		for (Face face : faces) {
			Edge[] faceEdges = face.getEdges();
			
			int[] prevEnds = faceEdges[faceEdges.length - 1].getEnds();
			int[] prevIndices = newEdges.get(prevEnds[0]).get(prevEnds[1]);
			for (Edge faceEdge : faceEdges) {
				int[] ends = faceEdge.getEnds();
				int[] edgeIndices = newEdges.get(ends[0]).get(ends[1]);
				
				Face pentagon = new Face(5);
				pentagon.setAllVertexPositions(centroidIndices.get(face), prevIndices[1],
						faceEdge.getEnds()[0], edgeIndices[0], edgeIndices[1]);
				gyroPolyhedron.addFace(pentagon);
				
				// Update previous indices for next iteration
				prevIndices = edgeIndices;
			}
		}
		
		gyroPolyhedron.setVertexNormalsToFaceNormals();
		return gyroPolyhedron;
	}
	
	/**
	 * Computes the "snub" polyhedron of this polyhedron. Essentially, each face
	 * is twisted, and each edge is replaced with two equilateral triangles. Each
	 * vertex also becomes a new face.
	 * 
	 * @return The snub polyhedron.
	 */
	public Polyhedron snub() {
		return this.gyro().dual();
	}
	
	/**
	 * Computes the "bevel" polyhedron of this polyhedron. Behaves similarly to
	 * expand in that there is one new face for each vertex and edge. Also known
	 * as cantitruncation.
	 * 
	 * @return The bevel polyhedron.
	 */
	public Polyhedron bevel() {
		return this.ambo().truncate();
	}
	
	/**
	 * Computes the "medial" polyhedron of this polyhedron. Adds vertices at the
	 * face centroids and edge midpoints. Each face is split into 2n triangles,
	 * where n is the number of vertices in the face. These triangles share a
	 * vertex at the face's centroid.
	 * 
	 * @return The medial polyhedron.
	 */
	public Polyhedron medial() {
		return this.bevel().dual();
	}
	
}
