package test;

import java.io.File;

import mesh.Mesh;
import mesh.polyhedra.archimedes.Cuboctahedron;
import mesh.polyhedra.archimedes.Icosidodecahedron;
import mesh.polyhedra.plato.Cube;
import mesh.polyhedra.plato.Dodecahedron;
import mesh.polyhedra.plato.Icosahedron;
import mesh.polyhedra.plato.Octahedron;
import mesh.polyhedra.plato.Tetrahedron;
import mesh.polyhedra.Polyhedron;

public class PolyhedraGenTest {

	public static void main(String[] args) {
		long start = System.currentTimeMillis();
//		Tetrahedron t = new Tetrahedron(1.0);
		Cube c = new Cube(1.0);
//		Octahedron o = new Octahedron(1.0);
//		Dodecahedron d = new Dodecahedron(1.0);
//		Icosahedron i = new Icosahedron(1.0);
//		Cuboctahedron co = new Cuboctahedron(1.0);
//		Icosidodecahedron id = new Icosidodecahedron(1.0);
//		Mesh test = c.truncate().canonicalize();
		Polyhedron test = c.snub().snub().snub().canonicalize(20);
		
		long end = System.currentTimeMillis();
		System.out.printf("Ran for %f seconds.", (end - start) / 1000.);
		File aout = new File("snubSnubSnubCube.obj");
		test.toObjFile(aout);
	}
	
}
