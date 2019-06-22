package mesh.polyhedra.plato;

import lombok.Getter;
import mesh.polyhedra.Polyhedron;

/**
 * Abstract class for a Platonic solid mesh centered at the origin.
 * 
 * @author Brian Yao
 */
public abstract class PlatonicSolid extends Polyhedron {

	@Getter
	private double edgeLength;
	
	/**
	 * @param edgeLength The length of each edge in this platonic solid.
	 */
	public PlatonicSolid(double edgeLength) {
		this.edgeLength = edgeLength;
	}

}
