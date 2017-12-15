package mesh.polyhedra;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Vector3d;

import mesh.Edge;
import mesh.Face;
import mesh.Mesh;
import mesh.struct.EdgeToAdjacentFace;
import mesh.struct.OrderedVertexToAdjacentEdge;
import mesh.struct.OrderedVertexToAdjacentFace;
import util.Canonicalize;

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
	 * Canonicalizes this polyhedron for the given number of iterations.
	 * See util.Canonicalize for more details.
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
	 * Canonicalizes this polyhedron until the change in position does not
	 * exceed the given threshold. That is, the algorithm terminates when no vertex
	 * moves more than the threshold after one iteration.
	 * 
	 * @param threshold The threshold for change in one iteration.
	 * @return The canonicalized version of this polyhedron.
	 */
	public Polyhedron canonicalize(double threshold) {
		Polyhedron canonicalized = this.clone();
		Canonicalize.adjust(canonicalized, threshold);
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
					newFace.setVertexIndex(j, newVertexIndex);
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
				newFace.setVertexIndex(i, newVertices.get(edge));
			}
			amboPolyhedron.addFace(newFace);
		}
		
		// Construct new faces in place of vertices of original
		OrderedVertexToAdjacentEdge ovtae = new OrderedVertexToAdjacentEdge(this);
		for (int i = 0 ; i < numVertexPositions() ; i++) {
			List<Edge> adjacentEdges = ovtae.getAdjacentEdges(i);
			Face newVertexFace = new Face(adjacentEdges.size());
			for (int j = 0 ; j < adjacentEdges.size() ; j++) {
				newVertexFace.setVertexIndex(j, newVertices.get(adjacentEdges.get(j)));
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
				int v = adjacentFaces[0].getVertexIndex(i);
				if (v == edgeEnds[0]) {
					int vnext = (i + 1) % adjacentFaces[0].numVertices();
					if (adjacentFaces[0].getVertexIndex(vnext) != edgeEnds[1]) {
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
			topFace.setAllVertexIndices(vface0, vface1, edgeEnds[1]);
			bottomFace.setAllVertexIndices(vface1, vface0, edgeEnds[0]);
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
				pentagon.setAllVertexIndices(centroidIndices.get(face), prevIndices[1],
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
	
	/**
	 * Computes the "exalt" polyhedron of this polyhedron. Equivalent to
	 * applying needle twice.
	 * 
	 * @return The exalt polyhedron.
	 */
	public Polyhedron exalt() {
		return this.needle().needle();
	}
	
	/**
	 * Computes the "yank" polyhedron of this polyhedron. Equivalent to
	 * applying zip twice.
	 * 
	 * @return The yank polyhedron.
	 */
	public Polyhedron yank() {
		return this.zip().zip();
	}
	
	/**
	 * Computes the "chamfer" polyhedron of this polyhedron. Truncates edges
	 * and replaces them with hexagonal faces.
	 * 
	 * @return The chamfer polyhedron.
	 */
	public Polyhedron chamfer() {
		return this.dual().subdivide().dual();
	}
	
	/**
	 * Computes the "subdivide" polyhedron of this polyhedron. Adds vertices at
	 * the midpoints of edges, and creates new triangular faces around original
	 * vertices. Equivalent to ambo without removing the original vertices.
	 * 
	 * @return The subdivide polyhedron.
	 */
	public Polyhedron subdivide() {
		Polyhedron subdividePolyhedron = new Polyhedron();
		for (Vector3d vertexPos : vertexPositions) {
			subdividePolyhedron.addVertexPosition(new Vector3d(vertexPos));
		}
		
		// Create new vertices, one at the midpoint of each edge
		Map<Edge, Integer> newVertices = new HashMap<>();
		int vertexIndex = subdividePolyhedron.numVertexPositions();
		for (Face face : faces) {
			// Create a new face for each face on this polyhedron
			Face newFace = new Face(face.numVertices());
			Edge[] edges = face.getEdges();
			for (int i = 0 ; i < face.numVertices() ; i++) {
				Edge edge = edges[i];
				if (newVertices.get(edge) == null) {
					Vector3d edgeMidpt = edge.midpoint();
					subdividePolyhedron.addVertexPosition(edgeMidpt);
					newVertices.put(edge, vertexIndex++);
				}
				newFace.setVertexIndex(i, newVertices.get(edge));
			}
			subdividePolyhedron.addFace(newFace);
		}
		
		// Create new faces for each vertex
		OrderedVertexToAdjacentEdge ovtae = new OrderedVertexToAdjacentEdge(this);
		for (int i = 0 ; i < this.numVertexPositions() ; i++) {
			List<Edge> adjacentEdges = ovtae.getAdjacentEdges(i);
			
			Edge prevEdge = adjacentEdges.get(adjacentEdges.size() - 1);
			for (Edge edge : adjacentEdges) {
				int prevVertex = newVertices.get(prevEdge);
				int currVertex = newVertices.get(edge);
				Face triangle = new Face(3);
				triangle.setAllVertexIndices(i, prevVertex, currVertex);
				
				subdividePolyhedron.addFace(triangle);
				
				// Update previous edge
				prevEdge = edge;
			}
		}
		
		subdividePolyhedron.setVertexNormalsToFaceNormals();
		return subdividePolyhedron;
	}
	
	/**
	 * Computes the "loft" polyhedron of this polyhedron. Adds a smaller
	 * version of this face, with n trapezoidal faces connecting the inner
	 * smaller version and the outer original version, where n is the number
	 * of vertices the face has.
	 * 
	 * @return The loft polyhedron.
	 */
	public Polyhedron loft() {
		Polyhedron loftPolyhedron = new Polyhedron();
		for (Vector3d vertexPos : vertexPositions) {
			loftPolyhedron.addVertexPosition(new Vector3d(vertexPos));
		}
		
		// Generate new vertices
		Map<Face, int[]> newVertices = new HashMap<>();
		int vertexIndex = loftPolyhedron.numVertexPositions();
		for (Face face : faces) {
			Face shrunk = new Face(face.numVertices());
			int[] newFaceVertices = new int[face.numVertices()];
			
			Vector3d centroid = face.centroid();
			for (int i = 0 ; i < face.numVertices() ; i++) {
				int index = face.getVertexIndex(i);
				Vector3d vertex = vertexPositions.get(index);
				Vector3d toCentroid = new Vector3d();
				toCentroid.sub(centroid, vertex);
				toCentroid.scale(0.3); // 0 < arbitrary scale factor < 1
				
				Vector3d newVertex = new Vector3d();
				newVertex.add(vertex, toCentroid);
				
				loftPolyhedron.addVertexPosition(newVertex);
				newFaceVertices[i] = vertexIndex;
				shrunk.setVertexIndex(i, vertexIndex);
				vertexIndex++;
			}
			
			newVertices.put(face, newFaceVertices);
			loftPolyhedron.addFace(shrunk);
		}
		
		// Generate new faces
		for (Face face : faces) {
			int[] newFaceVertices = newVertices.get(face);
			int prevIndex = face.getVertexIndex(face.numVertices() - 1);
			int newPrevIndex = newFaceVertices[face.numVertices() - 1];
			for (int i = 0 ; i < face.numVertices() ; i++) {
				int currIndex = face.getVertexIndex(i);
				int newCurrIndex = newFaceVertices[i];
				
				Face trapezoid = new Face(4);
				trapezoid.setAllVertexIndices(prevIndex, currIndex,
						newCurrIndex, newPrevIndex);
				loftPolyhedron.addFace(trapezoid);
				
				prevIndex = currIndex;
				newPrevIndex = newCurrIndex;
			}
		}
		
		loftPolyhedron.setVertexNormalsToFaceNormals();
		return loftPolyhedron;
	}
	
	/**
	 * Compute the "quinto" polyhedron of this polyhedron. Equivalent to an
	 * ortho but truncating the vertex at the center of original faces. This
	 * creates a small copy of the original face (but rotated).
	 * 
	 * @return The quinto polyhedron.
	 */
	public Polyhedron quinto() {
		Polyhedron quintoPolyhedron = new Polyhedron();
		for (Vector3d vertexPos : vertexPositions) {
			quintoPolyhedron.addVertexPosition(new Vector3d(vertexPos));
		}
		
		// Create new vertices at the midpoint of each edge and toward the
		// face's centroid
		Map<Integer, Map<Integer, int[]>> newVertices = new HashMap<>();
		int vertexIndex = quintoPolyhedron.numVertexPositions();
		for (Face face : faces) {
			Vector3d centroid = face.centroid();
			
			Edge[] edges = face.getEdges();
			for (int i = 0 ; i < edges.length ; i++) {
				int[] endsi = edges[i].getEnds();
				if (!newVertices.containsKey(endsi[0])) {
					newVertices.put(endsi[0], new HashMap<Integer, int[]>());
				}
				
				// The indices of the new vertices
				int[] newIndices = new int[2];
				
				Vector3d edgeMidpt = edges[i].midpoint();
				
				Vector3d fromCentroid = new Vector3d();
				fromCentroid.sub(edgeMidpt, centroid);
				fromCentroid.scale(0.3); // 0 < arbitrary scale factor < 1
				
				Vector3d newVertex = new Vector3d();
				newVertex.add(centroid, fromCentroid);
				
				// Check if the midpoint of the edge has already been added
				if (newVertices.containsKey(endsi[1]) &&
						newVertices.get(endsi[1]).containsKey(endsi[0])) {
					int midptVertex = newVertices.get(endsi[1]).get(endsi[0])[0];
					newIndices[0] = midptVertex;
				} else {
					quintoPolyhedron.addVertexPosition(edgeMidpt);
					newIndices[0] = vertexIndex++;
				}
				
				quintoPolyhedron.addVertexPosition(newVertex);
				newIndices[1] = vertexIndex++;
				newVertices.get(endsi[0]).put(endsi[1], newIndices);
			}
		}
		
		// Generate new faces
		for (Face face : faces) {
			Face centralFace = new Face(face.numVertices());
			Edge[] edges = face.getEdges();
			
			int[] prevEnds = edges[edges.length - 1].getEnds();
			int[] prevVertices = newVertices.get(prevEnds[0]).get(prevEnds[1]);
			int centralIndex = 0;
			for (Edge currEdge : edges) {
				int[] currEnds = currEdge.getEnds();
				int[] currVertices = newVertices.get(currEnds[0]).get(currEnds[1]);
				
				Face pentagon = new Face(5);
				pentagon.setAllVertexIndices(prevVertices[1], prevVertices[0],
						currEnds[0], currVertices[0], currVertices[1]);
				quintoPolyhedron.addFace(pentagon);
				
				centralFace.setVertexIndex(centralIndex++, currVertices[1]);
				
				// Update previous vertex indices
				prevVertices = currVertices;
			}
			quintoPolyhedron.addFace(centralFace);
		}
		
		quintoPolyhedron.setVertexNormalsToFaceNormals();
		return quintoPolyhedron;
	}
	
	/**
	 * Computes the "lace" polyhedron of this polyhedron. Like loft, but has
	 * on each face an antiprism of the original face instead of a prism.
	 * 
	 * @return The lace polyhedron.
	 */
	public Polyhedron lace() {
		Polyhedron lacePolyhedron = new Polyhedron();
		for (Vector3d vertexPos : vertexPositions) {
			lacePolyhedron.addVertexPosition(new Vector3d(vertexPos));
		}
		
		// Generate new vertices
		int vertexIndex = lacePolyhedron.numVertexPositions();
		for (Face face : faces) {
			Face twist = new Face(face.numVertices());
			int[] newFaceVertices = new int[face.numVertices()];
			
			Vector3d centroid = face.centroid();
			Edge[] edges = face.getEdges();
			for (int i = 0 ; i < edges.length ; i++) {
				Vector3d edgeMidpt = edges[i].midpoint();
				
				Vector3d fromCentroid = new Vector3d();
				fromCentroid.sub(edgeMidpt, centroid);
				fromCentroid.scale(0.3); // 0 < arbitrary scale factor < 1
				
				Vector3d newVertex = new Vector3d();
				newVertex.add(centroid, fromCentroid);
				
				lacePolyhedron.addVertexPosition(newVertex);
				newFaceVertices[i] = vertexIndex;
				twist.setVertexIndex(i, vertexIndex++);
			}
			
			lacePolyhedron.addFace(twist);
			
			// Generate triangle faces between twist and original
			for (int i = 0 ; i < edges.length ; i++) {
				int[] endsi = edges[i].getEnds();
				int currVertex = newFaceVertices[i];
				int nextVertex = newFaceVertices[(i + 1) % newFaceVertices.length];
				
				Face largeTriangle = new Face(3);
				Face smallTriangle = new Face(3);
				largeTriangle.setAllVertexIndices(currVertex, endsi[0], endsi[1]);
				smallTriangle.setAllVertexIndices(nextVertex, currVertex, endsi[1]);
				
				lacePolyhedron.addFaces(largeTriangle, smallTriangle);
			}
		}
		
		lacePolyhedron.setVertexNormalsToFaceNormals();
		return lacePolyhedron;
	}
	
	/**
	 * Computes the "stake" polyhedron of this polyhedron. Like lace, but
	 * instead of having a central face, there is a central vertex and 
	 * quadrilaterals around the center.
	 * 
	 * @return The stake polyhedron.
	 */
	public Polyhedron stake() {
		Polyhedron stakePolyhedron = new Polyhedron();
		for (Vector3d vertexPos : vertexPositions) {
			stakePolyhedron.addVertexPosition(new Vector3d(vertexPos));
		}
		
		// Generate new vertices
		int vertexIndex = stakePolyhedron.numVertexPositions();
		for (Face face : faces) {
			int[] newFaceVertices = new int[face.numVertices()];
			
			Vector3d centroid = face.centroid();
			stakePolyhedron.addVertexPosition(centroid);
			int centroidIndex = vertexIndex++;
			
			Edge[] edges = face.getEdges();
			for (int i = 0 ; i < edges.length ; i++) {
				Vector3d edgeMidpt = edges[i].midpoint();
				
				Vector3d fromCentroid = new Vector3d();
				fromCentroid.sub(edgeMidpt, centroid);
				fromCentroid.scale(0.3); // 0 < arbitrary scale factor < 1
				
				Vector3d newVertex = new Vector3d();
				newVertex.add(centroid, fromCentroid);
				
				stakePolyhedron.addVertexPosition(newVertex);
				newFaceVertices[i] = vertexIndex++;
			}
			
			// Generate the quads and triangles on this face
			for (int i = 0 ; i < edges.length ; i++) {
				int[] endsi = edges[i].getEnds();
				int currVertex = newFaceVertices[i];
				int nextVertex = newFaceVertices[(i + 1) % newFaceVertices.length];
				
				Face triangle = new Face(3);
				Face quad = new Face(4);
				triangle.setAllVertexIndices(currVertex, endsi[0], endsi[1]);
				quad.setAllVertexIndices(nextVertex, centroidIndex, currVertex, endsi[1]);
				
				stakePolyhedron.addFaces(triangle, quad);
			}
		}
		
		stakePolyhedron.setVertexNormalsToFaceNormals();
		return stakePolyhedron;
	}
	
}
