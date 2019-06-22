package mesh;

import javax.vecmath.Vector3d;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * A class for polygon meshes. The mesh is represented by an indexed structure,
 * in which each face's vertex positions are specified by indexes in a list of
 * all vertex positions used. The mesh is designed for the OBJ file format; any
 * mesh can be written to a file using toObjFile().
 * 
 * @author Brian Yao
 */
public class Mesh {

	protected List<Face> faces;
	
	protected List<Vector3d> vertexPositions;
	protected List<Vector3d> vertexNormals;
	
	/**
	 * Create an empty mesh.
	 */
	public Mesh() {
		faces = new ArrayList<>();
		vertexPositions = new ArrayList<>();
		vertexNormals = new ArrayList<>();
	}
	
	/**
	 * @return The number of vertex positions stored in this mesh.
	 */
	public int numVertexPositions() {
		return vertexPositions.size();
	}
	
	/**
	 * @return The set of all edges in this mesh.
	 */
	public Set<Edge> getEdges() {
		Set<Edge> edges = new HashSet<>();
		faces.forEach(face -> edges.addAll(Arrays.asList(face.getEdges())));
		return edges;
	}
	
	/**
	 * @return A list of all faces stored in this mesh.
	 */
	public List<Face> getFaces() {
		return faces;
	}
	
	/**
	 * Add the specified face to this mesh. This sets the face's mesh field to
	 * point to this mesh.
	 * 
	 * @param face The face to add to this mesh.
	 */
	public void addFace(Face face) {
		faces.add(face);
		face.setMesh(this);
	}
	
	/**
	 * Add the specified faces to this mesh. Calls addFace() on each face.
	 * 
	 * @param faces The faces to add to this mesh.
	 */
	public void addFaces(Face... faces) {
		Arrays.asList(faces).forEach(this::addFace);
	}
	
	/**
	 * Add the specified faces to this mesh. Calls addFace() on each face.
	 * 
	 * @param faces The faces to add to this mesh.
	 */
	public void addFaces(Collection<Face> faces) {
		faces.forEach(this::addFace);
	}
	
	/**
	 * Add a new vertex position (simply a point in 3D space). This is added
	 * to the end of the list of vertex positions, so it takes a new highest
	 * vertex position index.
	 * 
	 * @param position The 3D vertex position to store in this mesh's data.
	 */
	public void addVertexPosition(Vector3d position) {
		vertexPositions.add(position);
	}
	
	/**
	 * Add multiple vertex positions to this mesh. See addVertexPosition().
	 * 
	 * @param positions The array of vertex positions to add to this mesh.
	 */
	public void addVertexPositions(Vector3d ... positions) {
		vertexPositions.addAll(Arrays.asList(positions));
	}
	
	/**
	 * Add multiple vertex positions to this mesh. See addVertexPosition().
	 * 
	 * @param positions The collection of vertex positions to add to this mesh.
	 */
	public void addVertexPositions(Collection<Vector3d> positions) {
		vertexPositions.addAll(positions);
	}
	
	/**
	 * Set the vertex positions.
	 * 
	 * @param positions The new vertex positions.
	 */
	public void setVertexPositions(List<Vector3d> positions) {
		vertexPositions = positions;
	}
	
	/**
	 * Add a new vertex normal (simply a vector in 3D space). This is added
	 * to the end of the list of vertex normals, so it takes a new highest
	 * vertex normal index.
	 * 
	 * @param normal The 3D vertex normal to store in this mesh's data.
	 */
	public void addVertexNormal(Vector3d normal) {
		vertexNormals.add(normal);
	}
	
	/**
	 * Add multiple vertex normals to this mesh. See addVertexPosition().
	 * 
	 * @param normals The array of vertex normals to add to this mesh.
	 */
	public void addVertexNormals(Vector3d ... normals) {
		vertexNormals.addAll(Arrays.asList(normals));
	}
	
	/**
	 * Add multiple vertex normals to this mesh. See addVertexPosition().
	 * 
	 * @param normals The collection of vertex normals to add to this mesh.
	 */
	public void addVertexNormals(Collection<Vector3d> normals) {
		vertexNormals.addAll(normals);
	}
	
	/**
	 * Set the vertex normals.
	 * 
	 * @param normals The new vertex normals.
	 */
	public void setVertexNormals(List<Vector3d> normals) {
		vertexNormals = normals;
	}
	
	/**
	 * @return A list of all vertex positions stored in this mesh.
	 */
	public List<Vector3d> getVertexPositions() {
		return vertexPositions;
	}
	
	/**
	 * @return A list of all vertex normals stored in this mesh.
	 */
	public List<Vector3d> getVertexNormals() {
		return vertexNormals;
	}
	
	/**
	 * For each face of the mesh, set the normals at the face vertices to be
	 * the normal of the face itself (since the face is a polygon, it lies in
	 * a plane). These normals replace any previously stored vertex normals
	 * in the mesh.
	 */
	public void setVertexNormalsToFaceNormals() {
		vertexNormals.clear();
		for (int i = 0 ; i < faces.size() ; i++) {
			Vector3d faceNormal = faces.get(i).getFaceNormal();
			addVertexNormal(faceNormal);
			faces.get(i).setAllVertexIndicesTo(i);
		}
	}
	
	/**
	 * Compute the equivalent triangle mesh of this mesh. If this mesh is
	 * already a triangle mesh, this method returns a copy of it. Otherwise,
	 * in the new mesh, each polygon is divided into triangles using the
	 * divideToTriangle() method in Face.
	 * 
	 * @return A triangle mesh whose geometry is the same as this one's.
	 */
	public Mesh toTriangleMesh() {
		Mesh triangleMesh = new Mesh();
		triangleMesh.vertexPositions.addAll(vertexPositions);
		triangleMesh.vertexNormals.addAll(vertexNormals);
		for (Face face : faces) {
			if (face.numVertices() == 3) {
				triangleMesh.addFace(face);
			} else {
				Face[] triangles = face.divideIntoTriangles();
				triangleMesh.addFaces(triangles);
			}
		}
		return triangleMesh;
	}
	
	/**
	 * Write this mesh to a file on the OBJ format.
	 * 
	 * @param outFile The file to write to.
	 */
	public void toObjFile(File outFile) {
		try {
			outFile.createNewFile();
			FileWriter outWriter = new FileWriter(outFile);
			for (Vector3d vertexPos : vertexPositions) {
				outWriter.write(String.format("v %f %f %f\n", vertexPos.x, vertexPos.y, vertexPos.z));
			}
			for (Vector3d normal : vertexNormals) {
				outWriter.write(String.format("vn %f %f %f\n", normal.x, normal.y, normal.z));
			}
			for (Face face : faces) {
				outWriter.write(face.toOBJString());
				outWriter.write('\n');
			}
			outWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	@Override
	public String toString() {
		return String.format("Vertices:%s\nNormals:%s\nFaces:%s\n", vertexPositions, vertexNormals, faces);
	}
	
}
