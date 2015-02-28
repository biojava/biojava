/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 08.05.2004
 *
 */


package org.biojava.nbio.structure ;

import org.biojava.nbio.structure.jama.Matrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;



/** utility operations on Atoms, AminoAcids, etc.
 * <p>
 * Currently the
 * coordinates of an Atom are stored as an array of size 3
 * (double[3]). It would be more powerful to use Point3D from
 * javax.vecmath.  but unfortunately this is not a part of standard
 * java installations, since it comes with java3d . So to keep things
 * simple at the moment biojava does not depend on java3d.
 * </p>
 * @author Andreas Prlic
 * @since 1.4
 * @version %I% %G%
 */

public class Calc {

	private final static Logger logger = LoggerFactory.getLogger(Calc.class);

	// 180 / pi
	static final double RADIAN = 57.29577951 ;

	/** Radians per degree.
	 *
	 */
	public final static float radiansPerDegree = (float) (2 * Math.PI / 360);

	/** Degrees per radian.
	 *
	 */
	public final static float degreesPerRadian = (float) (360 / (2 * Math.PI));


	/**
	 * calculate distance between two atoms.
	 *
	 * @param a  an Atom object
	 * @param b  an Atom object
	 * @return a double	
	 */
	public static final double getDistance(Atom a, Atom b) {
		double x = a.getX() - b.getX();
		double y = a.getY() - b.getY();
		double z = a.getZ() - b.getZ();

		double s  = x * x  + y * y + z * z;

		return Math.sqrt(s);
	}


	/**
	 * Will calculate the *square* of distances between two atoms. This will be
	 * faster as it will not perform the final square root to get the actual
	 * distance. Use this if doing large numbers of distance comparisons - it is
	 * marginally faster than getDistance().
	 *
	 * @param a  an Atom object
	 * @param b  an Atom object
	 * @return a double
	 */
	public static double getDistanceFast(Atom a, Atom b) {
		double x = a.getX() - b.getX();
		double y = a.getY() - b.getY();
		double z = a.getZ() - b.getZ();

		return x * x  + y * y + z * z;
	}

	public static final Atom invert(Atom a) {
		double[] coords = new double[]{0.0,0.0,0.0} ;
		Atom zero = new AtomImpl();
		zero.setCoords(coords);
		return subtract(zero, a);
	}


	/** add two atoms ( a + b).
	 *
	 * @param a  an Atom object
	 * @param b  an Atom object
	 * @return an Atom object
	 */
	public static final Atom add(Atom a, Atom b){

		Atom c = new AtomImpl();
		c.setX( a.getX() + b.getX() );
		c.setY( a.getY() + b.getY() );
		c.setZ( a.getZ() + b.getZ() );

		return c ;
	}



	/** subtract two atoms ( a - b).
	 *
	 * @param a  an Atom object
	 * @param b  an Atom object
	 * @return n new Atom object representing the difference


	 */
	public static final Atom subtract(Atom a, Atom b) {
		Atom c = new AtomImpl();
		c.setX( a.getX() - b.getX() );
		c.setY( a.getY() - b.getY() );
		c.setZ( a.getZ() - b.getZ() );

		return c ;
	}

	/** Vector product (cross product).
	 *
	 * @param a  an Atom object
	 * @param b  an Atom object
	 * @return an Atom object
	 */
	public static final Atom vectorProduct(Atom a , Atom b){

		Atom c = new AtomImpl();
		c.setX( a.getY() * b.getZ() - a.getZ() * b.getY() ) ;
		c.setY( a.getZ() * b.getX() - a.getX() * b.getZ() ) ;
		c.setZ( a.getX() * b.getY() - a.getY() * b.getX() ) ;
		return c ;

	}

	/**
	 * Scalar product (dot product).
	 *
	 * @param a an Atom object
	 * @param b an Atom object
	 * @return a double
	 */
	public static final double scalarProduct(Atom a, Atom b) {
		return a.getX() * b.getX() + a.getY() * b.getY() + a.getZ() * b.getZ();
	}

	/** Gets the length of the vector (2-norm)
	 *
	 * @param a  an Atom object
	 * @return Square root of the sum of the squared elements
	 */
	public static final double amount(Atom a){
		return Math.sqrt(scalarProduct(a,a));
	}

	/** Get the angle between two vectors
	 *
	 * @param a  an Atom object
	 * @param b  an Atom object
	 * @return Angle between a and b, in degrees
	 */
	public static final double angle(Atom a, Atom b){

		double skalar;
		double angle;

		skalar = scalarProduct(a,b);

		angle = skalar/( amount(a) * amount (b) );
		angle = Math.acos(angle);
		angle = angle * RADIAN ;

		return angle;
	}

	/** return the unit vector of vector a .
	 *
	 * @param a  an Atom object
	 * @return an Atom object
	 */
	public static final Atom unitVector(Atom a) {
		double amount = amount(a) ;

		double[] coords = new double[3];

		coords[0] = a.getX() / amount ;
		coords[1] = a.getY() / amount ;
		coords[2] = a.getZ() / amount ;

		a.setCoords(coords);
		return a;

	}

	/**
	 * Calculate the torsion angle, i.e. the angle between the normal vectors of the 
	 * two plains a-b-c and b-c-d.
	 *
	 * @param a  an Atom object
	 * @param b  an Atom object
	 * @param c  an Atom object
	 * @param d  an Atom object
	 * @return a double
	 */
	public static final double torsionAngle(Atom a, Atom b, Atom c, Atom d) {

		Atom ab = subtract(a,b);
		Atom cb = subtract(c,b);
		Atom bc = subtract(b,c);
		Atom dc = subtract(d,c);

		Atom abc = vectorProduct(ab,cb);
		Atom bcd = vectorProduct(bc,dc);

		double angl = angle(abc,bcd) ;

		/* calc the sign: */
		Atom vecprod = vectorProduct(abc,bcd);
		double val = scalarProduct(cb,vecprod);
		if (val<0.0) angl = -angl ;

		return angl;
	}

	/**
	 * Calculate the phi angle.
	 *
	 * @param a  an AminoAcid object
	 * @param b  an AminoAcid object
	 * @return a double
	 * @throws StructureException if aminoacids not connected or if any of the 4 needed atoms missing
	 */
	public static final double getPhi(AminoAcid a, AminoAcid b) throws StructureException {

		if ( ! isConnected(a,b)){
			throw new StructureException("can not calc Phi - AminoAcids are not connected!") ;
		}

		Atom a_C  = a.getC();
		Atom b_N  = b.getN();
		Atom b_CA = b.getCA();
		Atom b_C  = b.getC();

		// C and N were checked in isConnected already
		if (b_CA==null) throw new StructureException("Can not calculate Phi, CA atom is missing");

		return torsionAngle(a_C,b_N,b_CA,b_C);
	}

	/**
	 * Calculate the psi angle.
	 *
	 * @param a  an AminoAcid object
	 * @param b  an AminoAcid object
	 * @return a double
	 * @throws StructureException if aminoacids not connected or if any of the 4 needed atoms missing
	 */
	public static final double getPsi(AminoAcid a, AminoAcid b) throws StructureException {
		if ( ! isConnected(a,b)) {
			throw new StructureException("can not calc Psi - AminoAcids are not connected!") ;
		}

		Atom a_N   = a.getN();
		Atom a_CA  = a.getCA();
		Atom a_C   = a.getC();
		Atom b_N   = b.getN();

		// C and N were checked in isConnected already
		if (a_CA==null) throw new StructureException("Can not calculate Psi, CA atom is missing");

		return torsionAngle(a_N,a_CA,a_C,b_N);

	}

	/**
	 * Test if two amino acids are connected, i.e.
	 * if the distance from C to N < 2.5 Angstrom.
	 *
	 * If one of the AminoAcids has an atom missing, returns false.
	 *
	 * @param a  an AminoAcid object
	 * @param b  an AminoAcid object
	 * @return true if ...
	 */
	public static final boolean isConnected(AminoAcid a, AminoAcid b) {
		Atom C = null ;
		Atom N = null;

		C = a.getC();
		N = b.getN();

		if ( C == null || N == null)
			return false;

		// one could also check if the CA atoms are < 4 A...
		double distance = getDistance(C,N);
		return distance < 2.5;
	}



	/** rotate a single atom aroud a rotation matrix.
	 * matrix must be a 3x3 matrix.
	 *
	 * If the matrix is indexed m[row][col], then the matrix will be
	 * pre-multiplied (y=atom*M)
	 * @param atom atom to be rotated
	 * @param m a rotation matrix represented as a double[3][3] array
	 */
	public static final void rotate(Atom atom, double[][] m){

		double x = atom.getX();
		double y = atom.getY() ;
		double z = atom.getZ();

		double nx = m[0][0] * x + m[0][1] * y +  m[0][2] * z ;
		double ny = m[1][0] * x + m[1][1] * y +  m[1][2] * z ;
		double nz = m[2][0] * x + m[2][1] * y +  m[2][2] * z ;


		atom.setX(nx);
		atom.setY(ny);
		atom.setZ(nz);
	}

	/** Rotate a structure.
	 *
	 * @param structure a Structure object
	 * @param rotationmatrix an array (3x3) of double representing the rotation matrix. 
	 * @throws StructureException ...
	 */
	public static final void rotate(Structure structure, double[][] rotationmatrix)
			throws StructureException
	{
		if ( rotationmatrix.length != 3 ) {
			throw new StructureException ("matrix does not have size 3x3 !");
		}
		AtomIterator iter = new AtomIterator(structure) ;
		while (iter.hasNext()) {
			Atom atom = iter.next() ;
			Calc.rotate(atom,rotationmatrix);
		}
	}

	/** rotate a structure .
	 *
	 * @param group a group object
	 * @param rotationmatrix an array (3x3) of double representing the rotation matrix. 
	 * @throws StructureException ...
	 */
	public static final void rotate(Group group, double[][] rotationmatrix) throws StructureException {

		if ( rotationmatrix.length != 3 ) {
			throw new StructureException ("matrix does not have size 3x3 !");
		}
		AtomIterator iter = new AtomIterator(group) ;
		while (iter.hasNext()) {
			Atom atom = null ;

			atom = iter.next() ;
			rotate(atom,rotationmatrix);

		}
	}

	/** Rotate an atom around a Matrix object.
	 *
	 * @param atom atom to be rotated
	 * @param m rotation matrix to be applied to the atom
	 */
	public static final void rotate(Atom atom, Matrix m){

		double x = atom.getX();
		double y = atom.getY() ;
		double z = atom.getZ();
		double[][] ad = new double[][]{{x,y,z}};

		Matrix am = new Matrix(ad);
		Matrix na = am.times(m);

		atom.setX(na.get(0,0));
		atom.setY(na.get(0,1));
		atom.setZ(na.get(0,2));

	}

	/** Rotate a group object.
	 *
	 * @param group  a group to be rotated
	 * @param m a Matrix object representing the rotation matrix
	 */
	public static final void rotate(Group group, Matrix m){

		AtomIterator iter = new AtomIterator(group) ;

		while (iter.hasNext()) {
			Atom atom = iter.next() ;
			rotate(atom,m);

		}

	}

	/** Rotate a structure object.
	 *
	 * @param structure the structure to be rotated
	 * @param m rotation matrix to be applied 
	 */
	public static final void rotate(Structure structure, Matrix m){

		AtomIterator iter = new AtomIterator(structure) ;

		while (iter.hasNext()) {
			Atom atom = iter.next() ;
			rotate(atom,m);

		}

	}

	/**
	 * Transforms an atom object, given a Matrix4d (i.e. the vecmath library 
	 * double-precision 4x4 rotation+translation matrix)
	 * @param atom
	 * @param m
	 */
	public static final void transform (Atom atom, Matrix4d m) {

		Point3d p = new Point3d(atom.getX(),atom.getY(),atom.getZ());
		m.transform(p);

		atom.setX(p.x);
		atom.setY(p.y);
		atom.setZ(p.z);
	}

	/**
	 * Transforms a group object, given a Matrix4d (i.e. the vecmath library 
	 * double-precision 4x4 rotation+translation matrix)
	 * @param group
	 * @param m
	 */
	public static final void transform (Group group, Matrix4d m) {
		AtomIterator iter = new AtomIterator(group) ;

		while (iter.hasNext()) {
			Atom atom = iter.next() ;
			transform(atom,m);

		}
	}

	/**
	 * Transforms a structure object, given a Matrix4d (i.e. the vecmath library 
	 * double-precision 4x4 rotation+translation matrix)
	 * @param structure
	 * @param m
	 */
	public static final void transform (Structure structure, Matrix4d m) {
		AtomIterator iter = new AtomIterator(structure) ;

		while (iter.hasNext()) {
			Atom atom = iter.next() ;
			transform(atom,m);

		}
	}

	/**
	 * Transforms a chain object, given a Matrix4d (i.e. the vecmath library 
	 * double-precision 4x4 rotation+translation matrix)
	 * @param chain
	 * @param m
	 */
	public static final void transform (Chain chain, Matrix4d m) {

		for (Group g:chain.getAtomGroups()) {
			for (Atom atom: g.getAtoms()) {
				transform(atom,m);
			}
		}
	}

	/**
	 * Translates an atom object, given a Vector3d (i.e. the vecmath library 
	 * double-precision 3-d vector)
	 * @param atom
	 * @param v
	 */
	public static final void translate (Atom atom, Vector3d v) {

		atom.setX(atom.getX()+v.x);
		atom.setY(atom.getY()+v.y);
		atom.setZ(atom.getZ()+v.z);
	}

	/**
	 * Translates a group object, given a Vector3d (i.e. the vecmath library 
	 * double-precision 3-d vector)
	 * @param group
	 * @param v
	 */
	public static final void translate (Group group, Vector3d v) {
		AtomIterator iter = new AtomIterator(group) ;

		while (iter.hasNext()) {
			Atom atom = iter.next() ;
			translate(atom,v);

		}
	}

	/**
	 * Translates a Structure object, given a Vector3d (i.e. the vecmath library 
	 * double-precision 3-d vector)
	 * @param structure
	 * @param v
	 */
	public static final void translate (Structure structure, Vector3d v) {
		AtomIterator iter = new AtomIterator(structure) ;

		while (iter.hasNext()) {
			Atom atom = iter.next() ;
			translate(atom,v);

		}
	}

	/** calculate structure + Matrix coodinates ... 
	 *
	 * @param s the structure to operate on
	 * @param matrix a Matrix object
	 */
	public static final void plus(Structure s, Matrix matrix){
		AtomIterator iter = new AtomIterator(s) ;
		Atom oldAtom = null;
		Atom rotOldAtom = null;
		while (iter.hasNext()) {
			Atom atom = null ;

			atom = iter.next() ;
			try {
				if ( oldAtom != null){
					logger.debug("before {}", getDistance(oldAtom,atom));
				}
			} catch (Exception e){
				logger.error("Exception: ", e);
			}
			oldAtom = (Atom)atom.clone();

			double x = atom.getX();
			double y = atom.getY() ;
			double z = atom.getZ();
			double[][] ad = new double[][]{{x,y,z}};

			Matrix am = new Matrix(ad);
			Matrix na = am.plus(matrix);

			double[] coords = new double[3] ;
			coords[0] = na.get(0,0);
			coords[1] = na.get(0,1);
			coords[2] = na.get(0,2);
			atom.setCoords(coords);
			try {
				if ( rotOldAtom != null){
					logger.debug("after {}", getDistance(rotOldAtom,atom));
				}
			} catch (Exception e){
				logger.error("Exception: ", e);
			}
			rotOldAtom  = (Atom) atom.clone();
		}

	}



	/** shift a structure with a vector.
	 *
	 * @param structure  a Structure object
	 * @param a          an Atom object representing a shift vector
	 */
	public static final void shift(Structure structure, Atom a ){

		AtomIterator iter = new AtomIterator(structure) ;
		while (iter.hasNext() ) {
			Atom atom = null ;

			atom = iter.next()  ;

			Atom natom = add(atom,a);
			double x = natom.getX();
			double y = natom.getY() ;
			double z = natom.getZ();
			atom.setX(x);
			atom.setY(y);
			atom.setZ(z);

		}
	}

	/** Shift a vector.
	 *
	 * @param a vector a
	 * @param b vector b
	 */
	public static final void shift(Atom a, Atom b){

		Atom natom = add(a,b);
		double x = natom.getX();
		double y = natom.getY() ;
		double z = natom.getZ();
		a.setX(x);
		a.setY(y);
		a.setZ(z);
	}

	/** Shift a Group with a vector.
	 *
	 * @param group   a group object
	 * @param a          an Atom object representing a shift vector
	 */
	public static final void shift(Group group , Atom a ){

		AtomIterator iter = new AtomIterator(group) ;
		while (iter.hasNext() ) {
			Atom atom = null ;

			atom = iter.next()  ;

			Atom natom = add(atom,a);
			double x = natom.getX();
			double y = natom.getY() ;
			double z = natom.getZ();
			atom.setX(x);
			atom.setY(y);
			atom.setZ(z);

		}
	}



	/** Returns the center  of mass of the set of atoms.
	 * @param atomSet a set of Atoms
	 * @return an Atom representing the Centroid of the set of atoms
	 */
	public static final Atom getCentroid(Atom[] atomSet){

		double[] coords = new double[3];

		coords[0] = 0;
		coords[1] = 0;
		coords[2] = 0 ;

		for (Atom a : atomSet) {
			coords[0] += a.getX();
			coords[1] += a.getY();
			coords[2] += a.getZ();
		}

		int n = atomSet.length;
		coords[0] = coords[0] / n;
		coords[1] = coords[1] / n;
		coords[2] = coords[2] / n;

		Atom vec = new AtomImpl();
		vec.setCoords(coords);
		return vec;

	}

	public static  Atom centerOfMass(Atom[] points) {
		Atom center = new AtomImpl();

		float totalMass = 0.0f;
		for (Atom a : points) {
			float mass = a.getElement().getAtomicMass();
			totalMass += mass;
			center = scaleAdd(mass, a, center);
		}

		center = scaleEquals(center, 1.0f/totalMass);
		return center;
	}

	/**
	 * Multiply elements of a by s (in place)
	 * @param a
	 * @param s
	 * @return the modified a
	 */
	public static Atom scaleEquals(Atom a, double s) {
		double x = a.getX();
		double y = a.getY();
		double z = a.getZ();

		x *= s;
		y *= s;
		z *= s;

		//Atom b = new AtomImpl();
		a.setX(x);
		a.setY(y);
		a.setZ(z);

		return a;
	}

	/**
	 * Multiply elements of a by s
	 * @param a
	 * @param s
	 * @return A new Atom with s*a
	 */
	public static Atom scale(Atom a, double s) {
		double x = a.getX();
		double y = a.getY();
		double z = a.getZ();

		Atom b = new AtomImpl();
		b.setX(x*s);
		b.setY(y*s);
		b.setZ(z*s);

		return b;
	}


	/**
	 * Perform linear transformation s*X+B, and store the result in b
	 * @param s Amount to scale x
	 * @param x Input coordinate
	 * @param b Vector to translate (will be modified)
	 * @return b, after modification
	 */
	public static Atom scaleAdd(double s, Atom x, Atom b){

		double xc = s*x.getX() + b.getX();
		double yc = s*x.getY() + b.getY();
		double zc = s*x.getZ() + b.getZ();

		//Atom a = new AtomImpl();
		b.setX(xc);
		b.setY(yc);
		b.setZ(zc);

		return b;
	}

	/** Returns the Vector that needs to be applied to shift a set of atoms
	 * to the Centroid.
	 * @param atomSet array of Atoms  
	 * @return the vector needed to shift the set of atoms to its geometric center
	 */
	public static final Atom getCenterVector(Atom[] atomSet){
		Atom centroid = getCentroid(atomSet);

		return getCenterVector(atomSet,centroid);

	}

	/** Returns the Vector that needs to be applied to shift a set of atoms
	 * to the Centroid, if the centroid is already known
	 * @param atomSet array of Atoms  
	 * @return the vector needed to shift the set of atoms to its geometric center
	 */
	public static final Atom getCenterVector(Atom[] atomSet, Atom centroid){


		double[] coords = new double[3];
		coords[0] = 0 - centroid.getX();
		coords[1] = 0 - centroid.getY();
		coords[2] = 0 - centroid.getZ();

		Atom shiftVec = new AtomImpl();
		shiftVec.setCoords(coords);
		return shiftVec;

	}


	/** Center the atoms at the Centroid. 
	 * @param atomSet a set of Atoms
	 * @return an Atom representing the Centroid of the set of atoms
	 * @throws StructureException
	 * */
	public static final Atom[] centerAtoms(Atom[] atomSet) throws StructureException {

		Atom centroid = getCentroid(atomSet);
		return centerAtoms(atomSet, centroid);
	}

	/** Center the atoms at the Centroid, if the centroid is already know.
	 * @param atomSet a set of Atoms
	 * @return an Atom representing the Centroid of the set of atoms
	 * @throws StructureException
	 * */
	public static final Atom[] centerAtoms(Atom[] atomSet, Atom centroid) throws StructureException {

		Atom shiftVector = getCenterVector(atomSet, centroid);

		Atom[] newAtoms = new AtomImpl[atomSet.length];

		for (int i =0 ; i < atomSet.length; i++){
			Atom a = atomSet[i];
			Atom n = add(a,shiftVector);
			newAtoms[i] = n ;
		}
		return newAtoms;
	}




	/** creates a virtual C-beta atom. this might be needed when working with GLY
	 *
	 * thanks to Peter Lackner for a python template of this method.
	 * @param amino the amino acid for which a "virtual" CB atom should be calculated 
	 * @return a "virtual" CB atom
	 * @throws StructureException
	 */
	public static final Atom createVirtualCBAtom(AminoAcid amino)
			throws StructureException{

		AminoAcid  ala = StandardAminoAcid.getAminoAcid("ALA");
		Atom aN  = ala.getN();
		Atom aCA = ala.getCA();
		Atom aC  = ala.getC();
		Atom aCB = ala.getCB();


		Atom[] arr1 = new Atom[3];
		arr1[0] = aN;
		arr1[1] = aCA;
		arr1[2] = aC;

		Atom[] arr2 = new Atom[3];
		arr2[0] = amino.getN();
		arr2[1] = amino.getCA();
		arr2[2] = amino.getC();

		// ok now we got the two arrays, do a SVD:

		SVDSuperimposer svd = new SVDSuperimposer(arr2,arr1);

		Matrix rotMatrix = svd.getRotation();
		Atom tranMatrix = svd.getTranslation();

		Calc.rotate(aCB,rotMatrix);

		Atom virtualCB = Calc.add(aCB,tranMatrix);
		virtualCB.setName("CB");

		return virtualCB;
	}


	/**Get euler angles for a matrix given in ZYZ convention.
	 * (as e.g. used by Jmol)
	 *
	 * @param m the rotation matrix
	 * @return the euler values for a rotation around Z, Y, Z in degrees...
	 */
	public static final double[] getZYZEuler(Matrix m) {
		double m22 = m.get(2,2);
		double rY = (float) Math.acos(m22) * degreesPerRadian;
		double rZ1, rZ2;
		if (m22 > .999d || m22 < -.999d) {
			rZ1 = Math.atan2(m.get(1,0),  m.get(1,1)) * degreesPerRadian;
			rZ2 = 0;
		} else {
			rZ1 = Math.atan2(m.get(2,1), -m.get(2,0)) * degreesPerRadian;
			rZ2 = Math.atan2(m.get(1,2),  m.get(0,2)) * degreesPerRadian;
		}
		return new double[] {rZ1,rY,rZ2};
	}


	/** Convert a rotation Matrix to Euler angles.
	 *   This conversion uses conventions as described on page:
	 *   http://www.euclideanspace.com/maths/geometry/rotations/euler/index.htm
	 *   Coordinate System: right hand
	 *   Positive angle: right hand
	 *   Order of euler angles: heading first, then attitude, then bank
	 *
	 * @param m the rotation matrix
	 * @return a array of three doubles containing the three euler angles in radians
	 */
	public static final double[] getXYZEuler(Matrix m){
		double heading, attitude, bank;

		// Assuming the angles are in radians.
		if (m.get(1,0) > 0.998) { // singularity at north pole
			heading = Math.atan2(m.get(0,2),m.get(2,2));
			attitude = Math.PI/2;
			bank = 0;

		} else if  (m.get(1,0) < -0.998) { // singularity at south pole
			heading = Math.atan2(m.get(0,2),m.get(2,2));
			attitude = -Math.PI/2;
			bank = 0;

		} else {
			heading = Math.atan2(-m.get(2,0),m.get(0,0));
			bank = Math.atan2(-m.get(1,2),m.get(1,1));
			attitude = Math.asin(m.get(1,0));
		}
		return new double[] { heading, attitude, bank };
	}



	/** This conversion uses NASA standard aeroplane conventions as described on page:
	 *   http://www.euclideanspace.com/maths/geometry/rotations/euler/index.htm
	 *   Coordinate System: right hand
	 *   Positive angle: right hand
	 *   Order of euler angles: heading first, then attitude, then bank.
	 *   matrix row column ordering:
	 *   [m00 m01 m02]
	 *   [m10 m11 m12]
	 *   [m20 m21 m22]
	 * @param heading in radians
	 * @param attitude  in radians
	 * @param bank  in radians
	 * @return the rotation matrix */
	public static final  Matrix matrixFromEuler(double heading, double attitude, double bank) {
		// Assuming the angles are in radians.
		double ch = Math.cos(heading);
		double sh = Math.sin(heading);
		double ca = Math.cos(attitude);
		double sa = Math.sin(attitude);
		double cb = Math.cos(bank);
		double sb = Math.sin(bank);

		Matrix m = new Matrix(3,3);
		m.set(0,0, ch * ca);
		m.set(0,1, sh*sb - ch*sa*cb);
		m.set(0,2, ch*sa*sb + sh*cb);
		m.set(1,0, sa);
		m.set(1,1, ca*cb);
		m.set(1,2, -ca*sb);
		m.set(2,0, -sh*ca);
		m.set(2,1, sh*sa*cb + ch*sb);
		m.set(2,2, -sh*sa*sb + ch*cb);

		return m;
	}


	/**
	 * Calculates the angle from centerPt to targetPt in degrees.
	 * The return should range from [0,360), rotating CLOCKWISE, 
	 * 0 and 360 degrees represents NORTH,
	 * 90 degrees represents EAST, etc...
	 *
	 * Assumes all points are in the same coordinate space.  If they are not, 
	 * you will need to call SwingUtilities.convertPointToScreen or equivalent 
	 * on all arguments before passing them  to this function.
	 *
	 * @param centerPt   Point we are rotating around.
	 * @param targetPt   Point we want to calculate the angle to.  
	 * @return angle in degrees.  This is the angle from centerPt to targetPt.
	 */
	public static double calcRotationAngleInDegrees(Atom centerPt, Atom targetPt)
	{
		// calculate the angle theta from the deltaY and deltaX values
		// (atan2 returns radians values from [-PI,PI])
		// 0 currently points EAST.  
		// NOTE: By preserving Y and X param order to atan2,  we are expecting 
		// a CLOCKWISE angle direction.  
		double theta = Math.atan2(targetPt.getY() - centerPt.getY(), targetPt.getX() - centerPt.getX());

		// rotate the theta angle clockwise by 90 degrees 
		// (this makes 0 point NORTH)
		// NOTE: adding to an angle rotates it clockwise.  
		// subtracting would rotate it counter-clockwise
		theta += Math.PI/2.0;

		// convert from radians to degrees
		// this will give you an angle from [0->270],[-180,0]
		double angle = Math.toDegrees(theta);

		// convert to positive range [0-360)
		// since we want to prevent negative angles, adjust them now.
		// we can assume that atan2 will not return a negative value
		// greater than one partial rotation
		if (angle < 0) {
			angle += 360;
		}

		return angle;
	}


	public static void main(String[] args){
		Atom a =new AtomImpl();
		a.setX(0);
		a.setY(0);
		a.setZ(0);

		Atom b = new AtomImpl();
		b.setX(1);
		b.setY(1);
		b.setZ(0);

		logger.info("Angle between atoms: ", calcRotationAngleInDegrees(a, b));
	}


	public static void rotate(Atom[] ca, Matrix matrix) {
		for (Atom atom : ca) Calc.rotate(atom, matrix);
	}
}


