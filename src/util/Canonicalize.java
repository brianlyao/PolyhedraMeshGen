package util;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Vector3d;

import mesh.Face;
import mesh.polyhedra.Polyhedron;

/**
 * A Java implementation of precisely the iterative algorithm for computing
 * canonical polyhedra designed by George W. Hart. All the code in this class
 * is based directly off of his work.
 * 
 * Further information on this algorithm and what it does is on Hart's website.
 * See this link: http://www.georgehart.com/canonical/canonical-supplement.html
 * 
 * @author Brian Yao
 */
public class Canonicalize {

	// Unused
	private static List<Vector3d> reciprocalVertices(Polyhedron poly) {
		List<Vector3d> newVertices = new ArrayList<>();
		
		List<Vector3d> vertexPositions = poly.getVertexPositions();
		for (Face face : poly.getFaces()) {
			// Initialize values which will be updated in the loop below
			Vector3d centroid = face.vertexAverage();
			Vector3d normalSum = new Vector3d();
			double avgEdgeDistance = 0.;
			
			// Retrieve the indices of the vertices defining this face
			int[] faceVertexIndices = face.getVertexIndices();
			
			// Keep track of the "previous" two vertices in CCW order
			int lastLastVertexIndex = faceVertexIndices[faceVertexIndices.length - 2];
			int lastVertexIndex = faceVertexIndices[faceVertexIndices.length - 1];
			for (int vertexIndex : faceVertexIndices) {
				Vector3d vertex = vertexPositions.get(vertexIndex);
				
				// Compute the normal of the plane defined by this vertex and
				// the previous two
				Vector3d lastlastVertex = vertexPositions.get(lastLastVertexIndex);
				Vector3d lastVertex = vertexPositions.get(lastVertexIndex);
				Vector3d v1 = new Vector3d(lastlastVertex);
				v1.sub(lastVertex);
				Vector3d v2 = new Vector3d(vertex);
				v2.sub(lastVertex);
				Vector3d normal = new Vector3d();
				normal.cross(v1, v2);
				normalSum.add(normal);
				
				// Compute distance from edge to origin
				avgEdgeDistance += Geometry.pointLineDist(new Vector3d(), lastlastVertex, lastVertex);
				
				// Update the previous vertices for the next iteration
				lastLastVertexIndex = lastVertexIndex;
				lastVertexIndex = vertexIndex;
			}
			
			normalSum.normalize();
			avgEdgeDistance /= faceVertexIndices.length;
			
			Vector3d resultingVector = new Vector3d();
			resultingVector.scale(centroid.dot(normalSum), normalSum);
			resultingVector.scale(1.0 / resultingVector.lengthSquared());
			resultingVector.scale((1.0 + avgEdgeDistance) / 2.0);
			newVertices.add(resultingVector);
		}
		
		return newVertices;
	}
	
	// Unused
	public static void canonicalize(Polyhedron poly, int numIterations) {
		Polyhedron dual = poly.dual();
		for (int i = 0 ; i < numIterations ; i++) {
			dual.setVertexPositions(reciprocalVertices(poly));
			poly.setVertexPositions(reciprocalVertices(dual));
		}
		poly.setVertexNormalsToFaceNormals();
	}
	
	/**
	 * Reflects the centers of faces across the unit sphere.
	 * 
	 * @param poly The polyhedron whose centers to invert.
	 * @return The list of inverted face centers.
	 */
	private static List<Vector3d> reciprocalCenters(Polyhedron poly) {
		List<Vector3d> faceCenters = new ArrayList<>();
		for (Face face : poly.getFaces()) {
			Vector3d newCenter = new Vector3d(face.vertexAverage());
			newCenter.scale(1.0 / newCenter.lengthSquared());
			faceCenters.add(newCenter);
		}
		return faceCenters;
	}
	
	/**
	 * Canonicalizes a polyhedron by adjusting its vertices iteratively.
	 * 
	 * @param poly          The polyhedron to canonicalize.
	 * @param numIterations The number of iterations to adjust for.
	 */
	public static void adjust(Polyhedron poly, int numIterations) {
		Polyhedron dual = poly.dual();
		for (int i = 0 ; i < numIterations ; i++) {
			List<Vector3d> newDualPositions = reciprocalCenters(poly);
			dual.setVertexPositions(newDualPositions);
			List<Vector3d> newPositions = reciprocalCenters(dual);
			poly.setVertexPositions(newPositions);
		}
		poly.setVertexNormalsToFaceNormals();
	}
	
}
