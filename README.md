# PolyhedraMeshGen

PolyhedraMeshGen is a Java 8+ application containing code for generating polyhedra meshes. The main focus will be on convex polyhedra with certain degrees of symmetry, but this could change in the future. In addition to basic polyhedra such as the Platonic and Archimedean solids, I plan to support the generation of related meshes such as prisms, antiprisms, and approximate spheres. Suggestions are welcome.

This is in part inspired by George W. Hart's work on generating polyhedra using Conway notation. I plan to include the same functionality here.

Currently, all meshes are in OBJ format since that is what I am familiar with. I may look into incorporating conversions to other formats.

# License

PolyhedraMeshGen is licensed under the GNU General Public License v3.0. See LICENSE.txt for details.

# External Dependencies

PolyhedraMeshGen relies on the following external library dependencies:

- [vecmath v1.5.2](https://mvnrepository.com/artifact/javax.vecmath/vecmath/1.5.2), a library containing implementations of vector and matrix math.

# Current Features

See the wiki for the features implemented in the most up-to-date version.