package math;

import javax.vecmath.Vector3d;

public class VectorMath {

	private static final double EPSILON = 1e-15;
	
	public static boolean isZero(Vector3d vector) {
		return Math.abs(vector.x) < EPSILON && Math.abs(vector.y) < EPSILON && Math.abs(vector.z) < EPSILON;
	}
	
}
