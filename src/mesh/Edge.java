package mesh;

import javax.vecmath.Vector3d;

/**
 * A class encapsulating a single edge between vertices of a mesh. The edge,
 * much like faces, is specified by the indexes of the vertex positions in
 * the Mesh. It requires exactly two vertices to define an edge by its
 * endpoint vertices.
 * 
 * The edge behaves like an unordered pair of integers. The equals() and
 * hashCode() methods have been overridden to ensure this behavior.
 * 
 * The edge stores a pointer to the mesh it belongs to, but this is null until
 * this edge is computed with some context. In particular, setMesh() needs to
 * be called for some methods to work properly. This is done implicitly for the
 * edges returned when calling getEdges() on a face.
 * 
 * @author Brian Yao
 */
public class Edge {

	private Mesh mesh;
	
	private int[] ends;
	
	/**
	 * Create an edge whose endpoints are the specified vertices.
	 * 
	 * @param vertex0 Endpoint vertex.
	 * @param vertex1 Endpoint vertex.
	 */
	public Edge(int vertex0, int vertex1) {
		ends = new int[2];
		ends[0] = vertex0;
		ends[1] = vertex1;
	}
	
	/**
	 * Compute the length of this edge. This requires this edge's mesh to be set.
	 * 
	 * @return The edge's length.
	 */
	public double length() {
		Vector3d diff = new Vector3d();
		diff.sub(mesh.vertexPositions.get(ends[0]), mesh.vertexPositions.get(ends[1]));
		return diff.length();
	}
	
	/**
	 * Compute the midpoint along this edge (arithmetic mean of the coordinates
	 * of the edge's endpoints). This requires this edge's mesh to be set.
	 * 
	 * @return The coordinate of the edge's midpoint.
	 */
	public Vector3d midpoint() {
		Vector3d midpt = new Vector3d();
		midpt.add(mesh.vertexPositions.get(ends[0]), mesh.vertexPositions.get(ends[1]));
		midpt.scale(0.5);
		return midpt;
	}
	
	/**
	 * Set the mesh that this edge points to. Just like with faces, this method
	 * needs to be called at some point for certain methods in this class and
	 * others to work properly, since Edge objects do not store any geometry;
	 * all geometry is stored in the Mesh.
	 * 
	 * @param mesh The mesh this edge will now point to.
	 */
	public void setMesh(Mesh mesh) {
		this.mesh = mesh;
	}
	
	@Override
	public boolean equals(Object obj) {
		Edge other = (Edge) obj;
		return (ends[0] == other.ends[0] && ends[1] == other.ends[1]) ||
			(ends[0] == other.ends[1] && ends[1] == other.ends[0]);
	}
	
	@Override
	public int hashCode() {
		return Integer.valueOf(ends[0]).hashCode() * Integer.valueOf(ends[1]).hashCode();
	}
	
}
