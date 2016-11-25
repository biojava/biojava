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
 */
package org.biojava.nbio.structure.xtal;

import javax.vecmath.*;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;

//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;

/**
 * A crystal cell's parameters.
 *
 * @author duarte_j
 *
 */
public class CrystalCell implements Serializable {

	private static final long serialVersionUID = 1L;
	//private static final Logger logger = LoggerFactory.getLogger(CrystalCell.class);

	public static final double MIN_VALID_CELL_SIZE = 10.0; // the minimum admitted for a crystal cell


	private double a;
	private double b;
	private double c;

	private double alpha;
	private double beta;
	private double gamma;

	private double alphaRad;
	private double betaRad;
	private double gammaRad;

	private double volume; // cached volume

	private double maxDimension; // cached max dimension

	private Matrix3d M; 	// cached basis change transformation matrix
	private Matrix3d Minv;  // cached basis change transformation matrix
	private Matrix3d Mtransp;
	private Matrix3d MtranspInv;

	public CrystalCell() {

	}

	public CrystalCell(double a, double b, double c, double alpha, double beta, double gamma){
		this.a = a;
		this.b = b;
		this.c = c;
		this.setAlpha(alpha); // initialises both alpha and alphaRad
		this.setBeta(beta);
		this.setGamma(gamma);
	}

	public double getA() {
		return a;
	}

	public void setA(double a) {
		this.a = a;
	}

	public double getB() {
		return b;
	}

	public void setB(double b) {
		this.b = b;
	}

	public double getC() {
		return c;
	}

	public void setC(double c) {
		this.c = c;
	}

	public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
		this.alphaRad = Math.toRadians(alpha);
	}

	public double getBeta() {
		return beta;
	}

	public void setBeta(double beta) {
		this.beta = beta;
		this.betaRad = Math.toRadians(beta);
	}

	public double getGamma() {
		return gamma;
	}

	public void setGamma(double gamma) {
		this.gamma = gamma;
		this.gammaRad = Math.toRadians(gamma);
	}

	/**
	 * Returns the volume of this unit cell.
	 * See http://en.wikipedia.org/wiki/Parallelepiped
	 * @return
	 */
	public double getVolume() {
		if (volume!=0) {
			return volume;
		}
		volume =  a*b*c*
		Math.sqrt(1-Math.cos(alphaRad)*Math.cos(alphaRad)-Math.cos(betaRad)*Math.cos(betaRad)-Math.cos(gammaRad)*Math.cos(gammaRad)
				+2.0*Math.cos(alphaRad)*Math.cos(betaRad)*Math.cos(gammaRad));

		return volume;
	}

	/**
	 * Get the index of a unit cell to which the query point belongs.
	 *
	 * <p>For instance, all points in the unit cell at the origin will return (0,0,0);
	 * Points in the unit cell one unit further along the `a` axis will return (1,0,0),
	 * etc.
	 * @param pt Input point (in orthonormal coordinates)
	 * @return A new point with the three indices of the cell containing pt
	 */
	public Point3i getCellIndices(Tuple3d pt) {
		Point3d p = new Point3d(pt);
		this.transfToCrystal(p);

		int x = (int)Math.floor(p.x);
		int y = (int)Math.floor(p.y);
		int z = (int)Math.floor(p.z);
		return new Point3i(x,y,z);
	}

	/**
	 * Converts the coordinates in pt so that they occur within the (0,0,0) unit cell
	 * @param pt
	 */
	public void transfToOriginCell(Tuple3d pt) {
		transfToCrystal(pt);

		// convert all coordinates to [0,1) interval
		pt.x = pt.x<0 ? (pt.x%1.0 + 1.0)%1.0 : pt.x%1.0;
		pt.y = pt.y<0 ? (pt.y%1.0 + 1.0)%1.0 : pt.y%1.0;
		pt.z = pt.z<0 ? (pt.z%1.0 + 1.0)%1.0 : pt.z%1.0;

		transfToOrthonormal(pt);
	}

	/**
	 * Converts a set of points so that the reference point falls in the unit cell.
	 *
	 * This is useful to transform a whole chain at once, allowing some of the
	 * atoms to be outside the unit cell, but forcing the centroid to be within it.
	 *
	 * @param points A set of points to transform (in orthonormal coordinates)
	 * @param reference The reference point, which is unmodified but which would
	 *    be in the unit cell were it to have been transformed. It is safe to
	 *    use a member of the points array here.
	 */
	public void transfToOriginCell(Tuple3d[] points, Tuple3d reference) {
		reference = new Point3d(reference);//clone
		transfToCrystal(reference);

		int x = (int)Math.floor(reference.x);
		int y = (int)Math.floor(reference.y);
		int z = (int)Math.floor(reference.z);

		for( Tuple3d point: points ) {
			transfToCrystal(point);
			point.x -= x;
			point.y -= y;
			point.z -= z;
			transfToOrthonormal(point);
		}
	}

	/**
	 *
	 * @param ops Set of operations in orthonormal coordinates
	 * @param reference Reference point, which should be in the unit cell after
	 *  each operation (also in orthonormal coordinates)
	 * @return A set of orthonormal operators with equivalent rotation to the
	 *  inputs, but with translation such that the reference point would fall
	 *  within the unit cell
	 */
	public Matrix4d[] transfToOriginCellOrthonormal(Matrix4d[] ops, Tuple3d reference) {
		// reference in crystal coords
		Point3d refXtal = new Point3d(reference);
		transfToCrystal(refXtal);

		// Convert to crystal coords
		Matrix4d[] opsXtal = new Matrix4d[ops.length];
		for(int i=0;i<ops.length;i++) {
			opsXtal[i] = transfToCrystal(ops[i]);
		}

		Matrix4d[] transformed = transfToOriginCellCrystal(opsXtal, refXtal);

		// Convert back to orthonormal
		for(int i=0;i<ops.length;i++) {
			transformed[i] = transfToOrthonormal(transformed[i]);
		}
		return transformed;
	}
	/**
	 *
	 * @param ops Set of operations in crystal coordinates
	 * @param reference Reference point, which should be in the unit cell after
	 *  each operation (also in crystal coordinates)
	 * @return A set of crystal operators with equivalent rotation to the
	 *  inputs, but with translation such that the reference point would fall
	 *  within the unit cell
	 */
	public Matrix4d[] transfToOriginCellCrystal(Matrix4d[] ops, Tuple3d reference) {
		Matrix4d[] transformed = new Matrix4d[ops.length];
		for(int j=0;j<ops.length;j++) {
			Matrix4d op = ops[j];

			// transform the reference point
			Point3d xXtal = new Point3d(reference);
			op.transform(xXtal);

			// Calculate unit cell of transformed reference
			int x = (int)Math.floor(xXtal.x);
			int y = (int)Math.floor(xXtal.y);
			int z = (int)Math.floor(xXtal.z);

			Matrix4d translation = new Matrix4d();
			translation.set(new Vector3d(-x,-y,-z));

			// Compose op with an additional translation operator
			translation.mul(op);

			Point3d ref2 = new Point3d(reference);
			translation.transform(ref2);

			transformed[j] = translation;
		}

		return transformed;
	}

	/**
	 * Transform given Matrix4d in crystal basis to the orthonormal basis using
	 * the PDB axes convention (NCODE=1)
	 * @param m
	 * @return
	 */
	public Matrix4d transfToOrthonormal(Matrix4d m) {
		Vector3d trans = new Vector3d(m.m03,m.m13,m.m23);
		transfToOrthonormal(trans);

		Matrix3d rot = new Matrix3d();
		m.getRotationScale(rot);
		// see Giacovazzo section 2.E, eq. 2.E.1
		// Rprime = MT-1 * R * MT
		rot.mul(this.getMTranspose());
		rot.mul(this.getMTransposeInv(),rot);

		return new Matrix4d(rot,trans,1.0);
	}

	/**
	 * Transforms the given crystal basis coordinates into orthonormal coordinates.
	 * e.g. transfToOrthonormal(new Point3d(1,1,1)) returns the orthonormal coordinates of the
	 * vertex of the unit cell.
	 * See Giacovazzo section 2.E, eq. 2.E.1 (or any linear algebra manual)
	 * @param v
	 */
	public void transfToOrthonormal(Tuple3d v) {
		// see Giacovazzo section 2.E, eq. 2.E.1
		getMTransposeInv().transform(v);
	}

	/**
	 * Transform given Matrix4d in orthonormal basis to the crystal basis using
	 * the PDB axes convention (NCODE=1)
	 * @param m
	 * @return
	 */
	public Matrix4d transfToCrystal(Matrix4d m) {
		Vector3d trans = new Vector3d(m.m03,m.m13,m.m23);
		transfToCrystal(trans);

		Matrix3d rot = new Matrix3d();
		m.getRotationScale(rot);
		// see Giacovazzo section 2.E, eq. 2.E.1 (the inverse equation)
		// R = MT * Rprime * MT-1
		rot.mul(this.getMTransposeInv());
		rot.mul(this.getMTranspose(),rot);

		return new Matrix4d(rot,trans,1.0);
	}

	/**
	 * Transforms the given orthonormal basis coordinates into crystal coordinates.
	 * See Giacovazzo eq 2.20 (or any linear algebra manual)
	 * @param v
	 */
	public void transfToCrystal(Tuple3d v) {
		getMTranspose().transform(v);
	}

	/**
	 * Returns the change of basis (crystal to orthonormal) transform matrix, that is
	 * M inverse in the notation of Giacovazzo.
	 * Using the PDB axes convention
	 * (CCP4 uses NCODE identifiers to distinguish the different conventions, the PDB one is called NCODE=1)
	 * The matrix is only calculated upon first call of this method, thereafter it is cached.
	 * See "Fundamentals of Crystallography" C. Giacovazzo, section 2.5 (eq 2.30)
	 *
	 * The non-standard orthogonalisation codes (NCODE for ccp4) are flagged in REMARK 285 after 2011's remediation
	 * with text: "THE ENTRY COORDINATES ARE NOT PRESENTED IN THE STANDARD CRYSTAL FRAME". There were only 148 PDB
	 * entries with non-standard code in 2011. See:
	 * http://www.wwpdb.org/documentation/2011remediation_overview-061711.pdf
	 * The SCALE1,2,3 records contain the correct transformation matrix (what Giacovazzo calls M matrix).
	 * In those cases if we calculate the M matrix following Giacovazzo's equations here, we get an entirely wrong one.
	 * Examples of PDB with non-standard orthogonalisation are 1bab and 1bbb.
	 * @return
	 */
	private Matrix3d getMInv() {
		if (Minv!=null) {
			return Minv;
		}

		// see eq. 2.30 Giacovazzo
		Minv =  new Matrix3d(                    this.a,                                           0,              0,
							  this.b*Math.cos(gammaRad),                   this.b*Math.sin(gammaRad),              0,
							  this.c*Math.cos(betaRad) , -this.c*Math.sin(betaRad)*getCosAlphaStar(),  1.0/getCstar());
		return Minv;
	}

	// another axes convention (from Giacovazzo) eq. 2.31b, not sure what would be the NCODE of this
//	private Matrix3d getMInv() {
//		if (Minv!=null) {
//			return Minv;
//		}
//		Minv = new Matrix3d( 1/getAstar(),  -(getCosGammaStar()/getSinGammaStar())/getAstar(), this.a*Math.cos(betaRad),
//				                        0,                   1/(getBstar()*getSinGammaStar()), this.b*Math.sin(alphaRad),
//				                        0,                                                  0,                  this.c);
//		return Minv;
//	}

	// and yet another axes convention (from Giacovazzo) eq. 2.31c, not sure what would be the NCODE of this
//	private Matrix3d getMInv() {
//		if (Minv!=null) {
//			return Minv;
//		}
//		Minv = new Matrix3d( alphaRad*Math.sin(gammaRad)*getSinBetaStar(), this.a*Math.cos(gammaRad), this.a*Math.sin(gammaRad)*getCosBetaStar(),
//				                                                        0,                    this.b,                                          0,
//				                                                        0, this.c*Math.cos(alphaRad),                   this.c*Math.sin(alphaRad));
//		return Minv;
//	}

	// relationships among direct and reciprocal lattice parameters
	// see Table 2.1 of chapter 2 of Giacovazzo
	@SuppressWarnings("unused")
	private double getAstar() {
		return (this.b*this.c*Math.sin(alphaRad))/getVolume();
	}

	@SuppressWarnings("unused")
	private double getBstar() {
		return (this.a*this.c*Math.sin(betaRad))/getVolume();
	}

	private double getCstar() {
		return (this.a*this.b*Math.sin(gammaRad))/getVolume();
	}

	private double getCosAlphaStar() {
		return (Math.cos(betaRad)*Math.cos(gammaRad)-Math.cos(alphaRad))/(Math.sin(betaRad)*Math.sin(gammaRad));
	}

	@SuppressWarnings("unused")
	private double getCosBetaStar() {
		return (Math.cos(alphaRad)*Math.cos(gammaRad)-Math.cos(betaRad))/(Math.sin(alphaRad)*Math.sin(gammaRad));
	}

	@SuppressWarnings("unused")
	private double getCosGammaStar() {
		return (Math.cos(alphaRad)*Math.cos(betaRad)-Math.cos(gammaRad))/(Math.sin(alphaRad)*Math.sin(betaRad));
	}

	@SuppressWarnings("unused")
	private double getSinAlphaStar() {
		return getVolume()/(this.a*this.b*this.c*Math.sin(betaRad)*Math.sin(gammaRad));
	}

	@SuppressWarnings("unused")
	private double getSinBetaStar() {
		return getVolume()/(this.a*this.b*this.c*Math.sin(alphaRad)*Math.sin(gammaRad));
	}

	@SuppressWarnings("unused")
	private double getSinGammaStar() {
		return getVolume()/(this.a*this.b*this.c*Math.sin(alphaRad)*Math.sin(betaRad));
	}

	/**
	 * Returns the change of basis (orthonormal to crystal) transform matrix, that is
	 * M in the notation of Giacovazzo.
	 * Using the PDB convention (NCODE=1).
	 * The matrix is only calculated upon first call of this method, thereafter it is cached.
	 * See "Fundamentals of Crystallography" C. Giacovazzo, section 2.5
	 * @return
	 */
	private Matrix3d getM() {
		if (M!=null){
			return M;
		}
		M = new Matrix3d();
		M.invert(getMInv());
		return M;
	}

	public Matrix3d getMTranspose() {
		if (Mtransp!=null){
			return Mtransp;
		}
		Matrix3d M = getM();
		Mtransp = new Matrix3d();
		Mtransp.transpose(M);
		return Mtransp;
	}

	private Matrix3d getMTransposeInv() {
		if (MtranspInv!=null){
			return MtranspInv;
		}
		MtranspInv = new Matrix3d();
		MtranspInv.invert(getMTranspose());
		return MtranspInv;
	}

	/**
	 * Gets the maximum dimension of the unit cell.
	 * @return
	 */
	public double getMaxDimension() {
		if (maxDimension!=0) {
			return maxDimension;
		}
		Point3d vert0 = new Point3d(0,0,0);
		Point3d vert1 = new Point3d(1,0,0);
		transfToOrthonormal(vert1);
		Point3d vert2 = new Point3d(0,1,0);
		transfToOrthonormal(vert2);
		Point3d vert3 = new Point3d(0,0,1);
		transfToOrthonormal(vert3);
		Point3d vert4 = new Point3d(1,1,0);
		transfToOrthonormal(vert4);
		Point3d vert5 = new Point3d(1,0,1);
		transfToOrthonormal(vert5);
		Point3d vert6 = new Point3d(0,1,1);
		transfToOrthonormal(vert6);
		Point3d vert7 = new Point3d(1,1,1);
		transfToOrthonormal(vert7);

		ArrayList<Double> vertDists = new ArrayList<Double>();
		vertDists.add(vert0.distance(vert7));
		vertDists.add(vert3.distance(vert4));
		vertDists.add(vert1.distance(vert6));
		vertDists.add(vert2.distance(vert5));
		maxDimension = Collections.max(vertDists);
		return maxDimension;
	}

	/**
	 * Given a scale matrix parsed from the PDB entry (SCALE1,2,3 records),
	 * checks that the matrix is a consistent scale matrix by comparing the
	 * cell volume to the inverse of the scale matrix determinant (tolerance of 1/100).
	 * If they don't match false is returned.
	 * See the PDB documentation for the SCALE record.
	 * See also last equation of section 2.5 of "Fundamentals of Crystallography" C. Giacovazzo
	 * @param scaleMatrix
	 * @return
	 */
	public boolean checkScaleMatrixConsistency(Matrix4d scaleMatrix) {

		double vol = getVolume();
		Matrix3d m = new Matrix3d();
		scaleMatrix.getRotationScale(m);

		// note we need to have a relaxed tolerance here as the PDB scale matrix is given with not such high precision
		// plus we don't want to have false positives, so we stay conservative
		double tolerance = vol/100.0;
		if ((Math.abs(vol - 1.0/m.determinant() )>tolerance)) {
			//System.err.println("Warning! SCALE matrix from PDB does not match 1/determinat == cell volume: "+
			//		String.format("vol=%6.3f  1/det=%6.3f",vol,1.0/m.determinant()));
			return false;
		}
		// this would be to check our own matrix, must always match!
		//if (!deltaComp(vol,1.0/getMTranspose().determinant())) {
		//	System.err.println("Our calculated SCALE matrix does not match 1/det=cell volume");
		//}

		return true;

	}

	/**
	 * Given a scale matrix parsed from a PDB entry (SCALE1,2,3 records),
	 * compares it to our calculated Mtranspose matrix to see if they coincide and
	 * returns true if they do.
	 * If they don't that means that the PDB entry is not in the standard
	 * orthogonalisation (NCODE=1 in ccp4).
	 * In 2011's remediation only 148 PDB entries were found not to be in
	 * a non-standard orthogonalisation. See:
	 * http://www.wwpdb.org/documentation/2011remediation_overview-061711.pdf
	 * For normal cases the scale matrix is diagonal without a translation component.
	 * Additionally the translation component of the SCALE matrix is also checked to
	 * make sure it is (0,0,0), if not false is return
	 * @param scaleMatrix
	 * @return
	 */
	public boolean checkScaleMatrix(Matrix4d scaleMatrix) {

		Matrix3d mtranspose = getMTranspose();
		for (int i=0;i<3;i++) {
			for (int j=0;j<3;j++) {
				if (!deltaComp(mtranspose.getElement(i, j),scaleMatrix.getElement(i, j))) {
					//System.out.println("Our value   ("+i+","+j+"): "+getM().getElement(i,j));
					//System.out.println("Their value ("+i+","+j+"): "+scaleMatrix.getElement(i,j));
					return false;
				}
			}
		}
		for (int i=0;i<3;i++) {
			if (!deltaComp(scaleMatrix.getElement(i, 3),0)) {
				return false;
			}
		}
		return true;
	}

	private boolean deltaComp(double d1, double d2) {
		if (Math.abs(d1-d2)<0.0001) return true;
		return false;
	}

	/**
	 * Checks whether the dimensions of this crystal cell are reasonable for protein
	 * crystallography: if all 3 dimensions are below {@value #MIN_VALID_CELL_SIZE} the cell
	 * is considered unrealistic and false returned
	 * @return
	 */
	public boolean isCellReasonable() {
		// this check is necessary mostly when reading PDB files that can contain the default 1 1 1 crystal cell
		// if we use that further things can go wrong, for instance for interface calculation
		// For instance programs like coot produce by default a 1 1 1 cell

		if (this.getA()<MIN_VALID_CELL_SIZE &&
				this.getB()<MIN_VALID_CELL_SIZE &&
				this.getC()<MIN_VALID_CELL_SIZE) {
			return false;
		}

		return true;

	}

	@Override
	public String toString() {
		return String.format("a%7.2f b%7.2f c%7.2f alpha%6.2f beta%6.2f gamma%6.2f", a, b, c, alpha, beta, gamma);
	}
}
