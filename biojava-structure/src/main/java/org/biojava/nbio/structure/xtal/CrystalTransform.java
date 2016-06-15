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

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;
import java.io.Serializable;

import static java.lang.Math.abs;


/**
 * Representation of a transformation in a crystal:
 * - a transformation id (each of the transformations in a space group, 0 to m)
 * - a crystal translation
 * The transformation matrix in crystal basis is stored, representing the basic
 * transformation together with the crystal translation.
 * Contains methods to check for equivalent transformations.
 *
 *
 * @author duarte_j
 *
 */
public class CrystalTransform implements Serializable {

	private static final long serialVersionUID = 1L;


	public static final Matrix4d IDENTITY = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);

	/**
	 * The space group to which this transform belongs
	 */
	private final SpaceGroup sg;

	/**
	 * The transform id corresponding to the SpaceGroup's transform indices.
	 * From 0 (identity) to m (m=number of symmetry operations of the space group)
	 * It is unique within the unit cell but equivalent units of different crystal unit cells
	 * will have same id
	 */
	private int transformId;

	/**
	 * The 4-dimensional matrix transformation in crystal basis.
	 * Note that the translational component of this matrix is not necessarily
	 * identical to crystalTranslation since some operators have fractional
	 * translations within the cell
	 */
	private Matrix4d matTransform;

	/**
	 * The crystal translation (always integer)
	 */
	private Point3i crystalTranslation;


	/**
	 * Creates a new CrystalTransform representing the identity transform
	 * in cell (0,0,0)
	 */
	public CrystalTransform(SpaceGroup sg) {
		this.sg = sg;
		this.transformId = 0;
		this.matTransform = (Matrix4d)IDENTITY.clone();
		this.crystalTranslation = new Point3i(0,0,0);
	}

	/**
	 * Represents the n-th transform
	 * @param sg
	 * @param transformId
	 */
	public CrystalTransform(SpaceGroup sg, int transformId) {
		this.sg = sg;
		this.transformId = transformId;
		if (sg==null && transformId==0) {
			this.matTransform = (Matrix4d)IDENTITY.clone();
		} else if (sg==null) {
			throw new IllegalArgumentException("Space Group cannot be null if transformId!=0");
		} else {
			this.matTransform = (Matrix4d)sg.getTransformation(transformId).clone();
		}
		this.crystalTranslation = new Point3i(0,0,0);
	}

	/**
	 * Copy constructor
	 * @param transform
	 */
	public CrystalTransform(CrystalTransform transform) {
		this.sg = transform.sg;
		this.transformId = transform.transformId;
		this.matTransform = new Matrix4d(transform.matTransform);
		this.crystalTranslation = new Point3i(transform.crystalTranslation);
	}

	public Matrix4d getMatTransform() {
		return matTransform;
	}

	public void setMatTransform(Matrix4d matTransform) {
		this.matTransform = matTransform;
	}

	public Point3i getCrystalTranslation() {
		return crystalTranslation;
	}

	public void translate(Point3i translation) {
		matTransform.m03 = matTransform.m03+translation.x;
		matTransform.m13 = matTransform.m13+translation.y;
		matTransform.m23 = matTransform.m23+translation.z;

		crystalTranslation.add(translation);

	}

	/**
	 * Returns true if the given CrystalTransform is equivalent to this one.
	 * Two crystal transforms are equivalent if one is the inverse of the other, i.e.
	 * their transformation matrices multiplication is equal to the identity.
	 * @param other
	 * @return
	 */
	public boolean isEquivalent(CrystalTransform other) {
		Matrix4d mul = new Matrix4d();
		mul.mul(this.matTransform,other.matTransform);

		if (mul.epsilonEquals(IDENTITY, 0.0001)) {
			return true;
		}
		return false;
	}

	/**
	 * Tells whether this transformation is a pure crystal lattice translation,
	 * i.e. no rotational component and an integer translation vector.
	 * @return
	 */
	public boolean isPureCrystalTranslation() {
		return (transformId==0 && (crystalTranslation.x!=0 || crystalTranslation.y!=0 || crystalTranslation.z!=0));
	}

	/**
	 * Tells whether this transformation is the identity: no rotation and no translation
	 * @return
	 */
	public boolean isIdentity() {
		return (transformId==0 && crystalTranslation.x==0 && crystalTranslation.y==0 && crystalTranslation.z==0);
	}

	/**
	 * Tells whether this transformation is a pure translation:
	 * either a pure crystal (lattice) translation or a fractional (within
	 * unit cell) translation: space groups Ixxx, Cxxx, Fxxx have operators
	 * with fractional translations within the unit cell.
	 * @return
	 */
	public boolean isPureTranslation() {
		if (isPureCrystalTranslation()) return true;
		if (SpaceGroup.deltaComp(matTransform.m00,1,SpaceGroup.DELTA) &&
			SpaceGroup.deltaComp(matTransform.m01,0,SpaceGroup.DELTA) &&
			SpaceGroup.deltaComp(matTransform.m02,0,SpaceGroup.DELTA) &&

			SpaceGroup.deltaComp(matTransform.m10,0,SpaceGroup.DELTA) &&
			SpaceGroup.deltaComp(matTransform.m11,1,SpaceGroup.DELTA) &&
			SpaceGroup.deltaComp(matTransform.m12,0,SpaceGroup.DELTA) &&

			SpaceGroup.deltaComp(matTransform.m20,0,SpaceGroup.DELTA) &&
			SpaceGroup.deltaComp(matTransform.m21,0,SpaceGroup.DELTA) &&
			SpaceGroup.deltaComp(matTransform.m22,1,SpaceGroup.DELTA) &&
			(	Math.abs(matTransform.m03-0.0)>SpaceGroup.DELTA ||
				Math.abs(matTransform.m13-0.0)>SpaceGroup.DELTA ||
				Math.abs(matTransform.m23-0.0)>SpaceGroup.DELTA)) {
			return true;
		}

		return false;
	}

	/**
	 * Tells whether this transformation contains a fractional translational
	 * component (whatever its rotational component). A fractional translation
	 * together with a rotation means a screw axis.
	 * @return
	 */
	public boolean isFractionalTranslation() {
		if ((Math.abs(matTransform.m03-crystalTranslation.x)>SpaceGroup.DELTA) ||
			(Math.abs(matTransform.m13-crystalTranslation.y)>SpaceGroup.DELTA) ||
			(Math.abs(matTransform.m23-crystalTranslation.z)>SpaceGroup.DELTA)) {
			return true;
		}
		return false;
	}

	/**
	 * Tells whether this transformation is a rotation disregarding the translational component,
	 * i.e. either pure rotation or screw rotation, but not improper rotation.
	 * @return
	 */
	public boolean isRotation() {
		// if no SG, that means a non-crystallographic entry (e.g. NMR). We return false
		if (sg==null) return false;

		// this also takes care of case <0 (improper rotations): won't be considered as rotations
		if (sg.getAxisFoldType(this.transformId)>1) return true;

		return false;
	}

	/**
	 * Returns the TransformType of this transformation: AU, crystal translation, fractional translation
	 * , 2 3 4 6-fold rotations, 2 3 4 6-fold screw rotations, -1 -3 -2 -4 -6 inversions/rotoinversions.
	 * @return
	 */
	public TransformType getTransformType() {

		// if no SG, that means a non-crystallographic entry (e.g. NMR). We return AU as type
		if (sg==null) return TransformType.AU;

		int foldType = sg.getAxisFoldType(this.transformId);
		boolean isScrewOrGlide = false;
		Vector3d translScrewComponent = getTranslScrewComponent();
		if (Math.abs(translScrewComponent.x-0.0)>SpaceGroup.DELTA ||
			Math.abs(translScrewComponent.y-0.0)>SpaceGroup.DELTA ||
			Math.abs(translScrewComponent.z-0.0)>SpaceGroup.DELTA) {

			isScrewOrGlide = true;
		}

		if (foldType>1) {

			if (isScrewOrGlide) {
				switch (foldType) {
				case 2:
					return TransformType.TWOFOLDSCREW;
				case 3:
					return TransformType.THREEFOLDSCREW;
				case 4:
					return TransformType.FOURFOLDSCREW;
				case 6:
					return TransformType.SIXFOLDSCREW;
				default:
					throw new NullPointerException("This transformation did not fall into any of the known types! This is most likely a bug.");
				}
			} else {
				switch (foldType) {
				case 2:
					return TransformType.TWOFOLD;
				case 3:
					return TransformType.THREEFOLD;
				case 4:
					return TransformType.FOURFOLD;
				case 6:
					return TransformType.SIXFOLD;
				default:
					throw new NullPointerException("This transformation did not fall into any of the known types! This is most likely a bug.");
				}
			}

		} else if (foldType<0) {
			switch (foldType) {
			case -1:
				return TransformType.ONEBAR;
			case -2:
				if (isScrewOrGlide) {
					return TransformType.GLIDE;
				}
				return TransformType.TWOBAR;
			case -3:
				return TransformType.THREEBAR;
			case -4:
				return TransformType.FOURBAR;
			case -6:
				return TransformType.SIXBAR;
			default:
				throw new NullPointerException("This transformation did not fall into any of the known types! This is most likely a bug.");
			}
		} else {
			if (isIdentity()) {
				return TransformType.AU;
			}
			if (isPureCrystalTranslation()) {
				return TransformType.XTALTRANSL;
			}
			if (isFractionalTranslation()) {
				return TransformType.CELLTRANSL;
			}
			throw new NullPointerException("This transformation did not fall into any of the known types! This is most likely a bug.");
		}

	}

	public Vector3d getTranslScrewComponent() {

		return getTranslScrewComponent(matTransform);

	}

	public int getTransformId() {
		return transformId;
	}

	public void setTransformId(int transformId) {
		this.transformId = transformId;
	}

	@Override
	public String toString() {
		return String.format("[%2d-(%s)]",transformId,toXYZString());
	}

	/**
	 * Expresses this transformation in terms of x,y,z fractional coordinates.
	 *
	 * Examples:
	 * @return
	 */
	public String toXYZString() {
		StringBuilder str = new StringBuilder();

		for(int i=0;i<3;i++) { //for each row
			boolean emptyRow = true;

			double coef; // TODO work with rational numbers


			// X
			coef = matTransform.getElement(i, 0);

			// Three cases for coef: zero, one, non-one
			if(abs(coef) > 1e-6 ) { // Non-zero
				if( abs( abs(coef)-1 ) < 1e-6 ) { // +/- 1
					if( coef < 0 ) {
						str.append("-");
					}
				} else {
					str.append(formatCoef(coef));
					str.append("*");
				}
				str.append("x");
				emptyRow = false;
			}

			// Y
			coef = matTransform.getElement(i, 1);

			if(abs(coef) > 1e-6 ) { // Non-zero
				if( abs( abs(coef)-1 ) < 1e-6 ) { // +/- 1
					if( coef < 0 ) {
						str.append("-");
					} else if( !emptyRow ) {
						str.append("+");
					}
				} else {
					if( !emptyRow && coef > 0) {
						str.append("+");
					}
					str.append(formatCoef(coef));
					str.append("*");
				}
				str.append("y");
				emptyRow = false;
			}

			// Z
			coef = matTransform.getElement(i, 2);

			if(abs(coef) > 1e-6 ) { // Non-zero
				if( abs( abs(coef)-1 ) < 1e-6 ) { // +/- 1
					if( coef < 0 ) {
						str.append("-");
					} else if( !emptyRow ) {
						str.append("+");
					}
				} else {
					if( !emptyRow && coef > 0) {
						str.append("+");
					}
					str.append(formatCoef(coef));
					str.append("*");
				}
				str.append("z");
				emptyRow = false;
			}

			// Intercept
			coef = matTransform.getElement(i, 3);

			if(abs(coef) > 1e-6 ) { // Non-zero
				if( !emptyRow && coef > 0) {
					str.append("+");
				}
				str.append(formatCoef(coef));
			}

			if(i<2) {
				str.append(",");
			}
		}


		return str.toString();
	}
	/**
	 * helper function to format simple fractions into rationals
	 * @param coef
	 * @return
	 */
	private String formatCoef(double coef) {
		double tol = 1e-6; // rounding tolerance

		// zero case
		if( Math.abs(coef) < tol) {
			return "0";
		}

		// integer case
		long num = Math.round(coef);
		if( Math.abs(num - coef) < tol) {
			return Long.toString(num);
		}

		// Other small cases
		for(int denom = 2; denom < 12; denom++ ) {
			num = Math.round(coef*denom);
			if( num - coef*denom < tol ) {
				return String.format("%d/%d",num, denom);
			}
		}

		// Give up and use floating point;
		return String.format("%.3f", coef);
	}

	/**
	 * Given a transformation matrix containing a rotation and translation returns the
	 * screw component of the rotation.
	 * See http://www.crystallography.fr/mathcryst/pdf/Gargnano/Aroyo_Gargnano_1.pdf
	 * @param m
	 * @return
	 */
	public static Vector3d getTranslScrewComponent(Matrix4d m) {

		int foldType = SpaceGroup.getRotAxisType(m);
		// For reference see:
		// http://www.crystallography.fr/mathcryst/pdf/Gargnano/Aroyo_Gargnano_1.pdf

		Vector3d transl = null;

		Matrix3d W =
				new Matrix3d(m.m00,m.m01,m.m02,
						m.m10,m.m11,m.m12,
						m.m20,m.m21,m.m22);

		if (foldType>=0) {

			// the Y matrix: Y = W^k-1 + W^k-2 ... + W + I  ; with k the fold type
			Matrix3d Y = new Matrix3d(1,0,0, 0,1,0, 0,0,1);
			Matrix3d Wk = new Matrix3d(1,0,0, 0,1,0, 0,0,1);

			for (int k=0;k<foldType;k++) {
				Wk.mul(W); // k=0 Wk=W, k=1 Wk=W^2, k=2 Wk=W^3, ... k=foldType-1, Wk=W^foldType
				if (k!=foldType-1) Y.add(Wk);
			}

			transl = new Vector3d(m.m03, m.m13, m.m23);
			Y.transform(transl);

			transl.scale(1.0/foldType);

		} else {

			if (foldType==-2) { // there are glide planes only in -2
				Matrix3d Y = new Matrix3d(1,0,0, 0,1,0, 0,0,1);
				Y.add(W);

				transl = new Vector3d(m.m03, m.m13, m.m23);
				Y.transform(transl);

				transl.scale(1.0/2.0);
			} else { // for -1, -3, -4 and -6 there's nothing to do: fill with 0s
				transl = new Vector3d(0,0,0);
			}
		}

		return transl;
	}
}
