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
package org.biojava.nbio.structure.align.util;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.structure.geometry.Matrices;
import org.biojava.nbio.structure.jama.Matrix;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import java.io.StringWriter;

/**
 * Calculates the rotation axis for an alignment
 *
 * <p>A superposition of two structures is generally represented as a rotation
 * matrix plus a translation vector. However, it can also be represented as an
 * axis of rotation plus some translation.
 *
 * <p>This class calculates the rotation axis and stores it as four properties:
 * <ul><li>A unit vector parallel to the rotation axis ({@link #getRotationAxis()})
 * <li>The angle of rotation ({@link #getAngle()})
 * <li>A point on the rotation axis ({@link #getRotationPos()})
 * <li>Some translation parallel to the axis ({@link #getScrewTranslation()})
 * </ul>
 *
 * <p>The axis of rotation is poorly defined and numerically unstable for small
 * angles. Therefore it's direction is left as null for angles less than
 * {@link #MIN_ANGLE}.
 *
 * @author Spencer Bliven
 *
 */
public final class RotationAxis {

	/**
	 * Minimum angle to calculate rotation axes. 5 degrees.
	 */
	static final double MIN_ANGLE = Math.toRadians(5.);

	private double theta;
	private Atom rotationAxis; // axis of rotation
	private Atom rotationPos; // a point on the axis of rotation
	private Atom screwTranslation; //translation parallel to the axis of rotation
	private Atom otherTranslation; // translation perpendicular to the axis of rotation

	/**
	 * The rotation angle
	 * @return the angle, in radians
	 */
	public double getAngle() {
		return theta;
	}

	/**
	 * Get a unit vector along the rotation axis
	 * @return rotationAxis
	 */
	public Atom getRotationAxis() {
		return rotationAxis;
	}

	/**
	 * Returns the rotation axis and angle in a single javax.vecmath.AxisAngle4d object
	 * @return
	 */
	public AxisAngle4d getAxisAngle4d() {
		return new AxisAngle4d(rotationAxis.getX(),rotationAxis.getY(),rotationAxis.getZ(),theta);
	}

	/**
	 * Get a position on the rotation axis.
	 *
	 * Specifically, project the origin onto the rotation axis
	 * @return rotationPos
	 */
	public Atom getRotationPos() {
		return rotationPos;
	}

	/**
	 * Get the component of translation parallel to the axis of rotation
	 * @return
	 */
	public Atom getScrewTranslation() {
		return screwTranslation;
	}

	public Vector3d getVector3dScrewTranslation() {
		return new Vector3d(screwTranslation.getX(),screwTranslation.getY(),screwTranslation.getZ());
	}
	
	public double getTranslation() {
		return Calc.amount(screwTranslation);
	}

	/**
	 * Calculate the rotation axis for the first block of an AFPChain
	 * @param afpChain
	 * @throws StructureException
	 * @throws NullPointerException if afpChain does not contain a valid rotation matrix and shift vector
	 */
	public RotationAxis(AFPChain afpChain) throws StructureException {
		if(afpChain.getAlnLength() < 1) {
			throw new StructureException("No aligned residues");
		}
		init(afpChain.getBlockRotationMatrix()[0],afpChain.getBlockShiftVector()[0]);
	}

	/**
	 * Create a rotation axis from a vector, a point, and an angle.
	 *
	 * The result will be a pure rotation, with no screw component.
	 * @param axis A vector parallel to the axis of rotation
	 * @param pos A point on the axis of rotation
	 * @param theta The angle to rotate (radians)
	 */
	public RotationAxis(Atom axis, Atom pos, double theta) {
		this.rotationAxis = Calc.unitVector(axis);
		this.rotationPos = (Atom) pos.clone();
		this.theta = theta;
		this.screwTranslation = new AtomImpl(); //zero
		this.otherTranslation = null; //deprecated
	}

	/**
	 * Determine the location of the rotation axis based on a rotation matrix and a translation vector
	 * @param rotation
	 * @param translation
	 */
	public RotationAxis(Matrix rotation, Atom translation) {
		init(rotation, translation);
	}

	/**
	 * Create a rotation axis from a Matrix4d containing a rotational
	 * component and a translational component.
	 *
	 * @param transform
	 */
	public RotationAxis(Matrix4d transform) {

		Matrix rot = Matrices.getRotationJAMA(transform);
		Atom transl = Calc.getTranslationVector(transform);
		init(rot,transl);
	}

	/**
	 * Get the rotation matrix corresponding to this axis
	 * @return A 3x3 rotation matrix
	 */
	public Matrix getRotationMatrix() {
		return getRotationMatrix(theta);
	}

	/**
	 * Get the rotation matrix corresponding to a rotation about this axis
	 * @param theta The amount to rotate
	 * @return A 3x3 rotation matrix
	 */
	public Matrix getRotationMatrix(double theta) {
		if( rotationAxis == null) {
			// special case for pure translational axes
			return Matrix.identity(3, 3);
		}
		double x = rotationAxis.getX();
		double y = rotationAxis.getY();
		double z = rotationAxis.getZ();
		double cos = Math.cos(theta);
		double sin = Math.sin(theta);
		double com = 1 - cos;
		return new Matrix(new double[][] {
				{com*x*x + cos, com*x*y+sin*z, com*x*z+-sin*y},
				{com*x*y-sin*z, com*y*y+cos, com*y*z+sin*x},
				{com*x*z+sin*y, com*y*z-sin*x, com*z*z+cos},
		});
	}

	/**
	 * Returns the rotation order o that gives the lowest value of {@code |2PI / o - theta},
	 * given that the value is strictly lower than {@code threshold}, for orders {@code o=1,...,maxOrder}.
	 */
	public int guessOrderFromAngle(double threshold, int maxOrder) {
		double bestDelta = threshold;
		int bestOrder = 1;
		for (int order = 2; order < maxOrder; order++) {
			double delta = Math.abs(2 * Math.PI / order - theta);
			if (delta < bestDelta) {
				bestOrder = order;
				bestDelta = delta;
			}
		}
		return bestOrder;
	}


	/**
	 * Returns a matrix that describes both rotation and translation.
	 */
	public Matrix getFullMatrix() {
		return null; // TODO, easy
	}

	/**
	 * Initialize variables
	 *
	 * @param rotation
	 * @param translation
	 */
	private void init(Matrix rotation, Atom translation) {
		if(rotation.getColumnDimension() != 3 || rotation.getRowDimension() != 3) {
			throw new IllegalArgumentException("Expected 3x3 rotation matrix");
		}


		// Calculate angle
		double c = (rotation.trace()-1)/2.0; //=cos(theta)
		// c is sometimes slightly out of the [-1,1] range due to numerical instabilities
		if( -1-1e-8 < c && c < -1 ) c = -1;
		if( 1+1e-8 > c && c > 1 ) c = 1;
		if( -1 > c || c > 1 ) {
			throw new IllegalArgumentException("Input matrix is not a valid rotation matrix.");
		}
		this.theta = Math.acos(c);

		if(theta < MIN_ANGLE) {
			calculateTranslationalAxis(rotation,translation);
		} else {
			calculateRotationalAxis(rotation, translation, c);
		}
	}

	/**
	 * Calculate the rotation axis for the normal case, where there is a
	 * significant rotation angle
	 * @param rotation
	 * @param translation
	 * @param c
	 */
	private void calculateRotationalAxis(Matrix rotation, Atom translation,
			double c) {
		// Calculate magnitude of rotationAxis components, but not signs
		double sum=0;
		double[] rotAx = new double[3];
		for(int i=0;i<3;i++) {
			rotAx[i] = Math.sqrt(rotation.get(i, i)-c);
			sum+=rotAx[i]*rotAx[i];
		}
		for(int i=0;i<3;i++) {
			rotAx[i] /= Math.sqrt(sum);
		}

		// Now determine signs
		double d0 = rotation.get(2,1)-rotation.get(1,2); //=2u[0]*sin(theta)
		double d1 = rotation.get(0,2)-rotation.get(2,0); //=2u[1]*sin(theta)
		double d2 = rotation.get(1,0)-rotation.get(0,1); //=2u[2]*sin(theta)

		double s12 = rotation.get(2,1)+rotation.get(1,2); //=2*u[1]*u[2]*(1-cos(theta))
		double s02 = rotation.get(0,2)+rotation.get(2,0); //=2*u[0]*u[2]*(1-cos(theta))
		double s01 = rotation.get(1,0)+rotation.get(0,1); //=2*u[0]*u[1]*(1-cos(theta))

		//Only use biggest d for the sign directly, for numerical stability around 180deg
		if( Math.abs(d0) < Math.abs(d1) ) { // not d0
			if( Math.abs(d1) < Math.abs(d2) ) { //d2
				if(d2>=0){ //u[2] positive
					if( s02 < 0 ) rotAx[0] = -rotAx[0];
					if( s12 < 0 ) rotAx[1] = -rotAx[1];
				} else { //u[2] negative
					rotAx[2] = -rotAx[2];
					if( s02 >= 0 ) rotAx[0] = -rotAx[0];
					if( s12 >= 0 ) rotAx[1] = -rotAx[1];
				}
			} else { //d1
				if(d1>=0) {//u[1] positive
					if( s01 < 0) rotAx[0] = -rotAx[0];
					if( s12 < 0) rotAx[2] = -rotAx[2];
				} else { //u[1] negative
					rotAx[1] = -rotAx[1];
					if( s01 >= 0) rotAx[0] = -rotAx[0];
					if( s12 >= 0) rotAx[2] = -rotAx[2];
				}
			}
		} else { // not d1
			if( Math.abs(d0) < Math.abs(d2) ) { //d2
				if(d2>=0){ //u[2] positive
					if( s02 < 0 ) rotAx[0] = -rotAx[0];
					if( s12 < 0 ) rotAx[1] = -rotAx[1];
				} else { //u[2] negative
					rotAx[2] = -rotAx[2];
					if( s02 >= 0 ) rotAx[0] = -rotAx[0];
					if( s12 >= 0 ) rotAx[1] = -rotAx[1];
				}
			} else { //d0
				if(d0>=0) { //u[0] positive
					if( s01 < 0 ) rotAx[1] = -rotAx[1];
					if( s02 < 0 ) rotAx[2] = -rotAx[2];
				} else { //u[0] negative
					rotAx[0] = -rotAx[0];
					if( s01 >= 0 ) rotAx[1] = -rotAx[1];
					if( s02 >= 0 ) rotAx[2] = -rotAx[2];
				}
			}
		}

		rotationAxis = new AtomImpl();
		rotationAxis.setCoords(rotAx);

		// Calculate screw = (rotationAxis dot translation)*u
		double dotProduct = Calc.scalarProduct(rotationAxis, translation);

		screwTranslation = Calc.scale(rotationAxis, dotProduct);
		otherTranslation = Calc.subtract(translation, screwTranslation);

		Atom hypot = Calc.vectorProduct(otherTranslation,rotationAxis);
		Calc.scaleEquals(hypot,.5/Math.tan(theta/2.0));

		// Calculate rotation axis position
		rotationPos = Calc.scaleAdd(.5,otherTranslation, hypot);
	}

	/**
	 * Handle cases with small angles of rotation
	 * @param rotation
	 * @param translation
	 */
	private void calculateTranslationalAxis(Matrix rotation, Atom translation) {
		// set axis parallel to translation
		rotationAxis = Calc.scale(translation, 1./Calc.amount(translation));

		// position is undefined
		rotationPos = null;

		screwTranslation = translation;
		otherTranslation = new AtomImpl();
		otherTranslation.setCoords(new double[] {0,0,0});
	}

	/**
	 * Returns a Jmol script which will display the axis of rotation. This
	 * consists of a cyan arrow along the axis, plus an arc showing the angle
	 * of rotation.
	 * <p>
	 * As the rotation angle gets smaller, the axis of rotation becomes poorly
	 * defined and would need to get farther and farther away from the protein.
	 * This is not particularly useful, so we arbitrarily draw it parallel to
	 * the translation and omit the arc.
	 * @param atoms Some atoms from the protein, used for determining the bounds
	 *  	  of the axis.
	 *
	 * @return The Jmol script, suitable for calls to
	 * {@link org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol#evalString() jmol.evalString()}
	 */
	public String getJmolScript(Atom[] atoms){
		return getJmolScript(atoms, 0);
	}

	/**
	 * Find a segment of the axis that covers the specified set of atoms.
	 * <p>
	 * Projects the input atoms onto the rotation axis and returns the bounding
	 * points.
	 * <p>
	 * In the case of a pure translational axis, the axis location is undefined
	 * so the center of mass will be used instead.
	 * @param atoms
	 * @return two points defining the axis segment
	 */
	public Pair<Atom> getAxisEnds(Atom[] atoms) {
		// Project each Atom onto the rotation axis to determine limits
		double min, max;
		min = max = Calc.scalarProduct(rotationAxis,atoms[0]);
		for(int i=1;i<atoms.length;i++) {
			double prod = Calc.scalarProduct(rotationAxis,atoms[i]);
			if(prod<min) min = prod;
			if(prod>max) max = prod;
		}
		double uLen = Calc.scalarProduct(rotationAxis,rotationAxis);// Should be 1, but double check
		min/=uLen;
		max/=uLen;
		
		// Project the origin onto the axis. If the axis is undefined, use the center of mass
		Atom axialPt;
		if(rotationPos == null) {
			Atom center = Calc.centerOfMass(atoms);

			// Project center onto the axis
			Atom centerOnAxis = Calc.scale(rotationAxis, Calc.scalarProduct(center, rotationAxis));

			// Remainder is projection of origin onto axis
			axialPt = Calc.subtract(center, centerOnAxis);

		} else {
			axialPt = rotationPos;
		}

		// Find end points of the rotation axis to display
		Atom axisMin = (Atom) axialPt.clone();
		Calc.scaleAdd(min, rotationAxis, axisMin);
		Atom axisMax = (Atom) axialPt.clone();
		Calc.scaleAdd(max, rotationAxis, axisMax);

		return new Pair<>(axisMin, axisMax);
	}
	/**
	 * Returns a Jmol script which will display the axis of rotation. This
	 * consists of a cyan arrow along the axis, plus an arc showing the angle
	 * of rotation.
	 * <p>
	 * As the rotation angle gets smaller, the axis of rotation becomes poorly
	 * defined and would need to get farther and farther away from the protein.
	 * This is not particularly useful, so we arbitrarily draw it parallel to
	 * the translation and omit the arc.
	 * @param atoms Some atoms from the protein, used for determining the bounds
	 *  	  of the axis.
	 * @param axisID in case of representing more than one axis in the same jmol
	 * 		  panel, indicate the ID number.
	 *
	 * @return The Jmol script, suitable for calls to
	 * {@link org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol#evalString() jmol.evalString()}
	 */
	public String getJmolScript(Atom[] atoms, int axisID){
		final double width=.5;// width of JMol object
		final String axisColor = "yellow"; //axis color
		final String screwColor = "orange"; //screw translation color

		Pair<Atom> endPoints = getAxisEnds(atoms);
		Atom axisMin = endPoints.getFirst();
		Atom axisMax = endPoints.getSecond();
		
		StringWriter result = new StringWriter();

		// set arrow heads to a reasonable length
		result.append("set defaultDrawArrowScale 2.0;");

		// draw axis of rotation
		result.append(
				String.format("draw ID rot"+axisID+" CYLINDER {%f,%f,%f} {%f,%f,%f} WIDTH %f COLOR %s ;",
						axisMin.getX(),axisMin.getY(),axisMin.getZ(),
						axisMax.getX(),axisMax.getY(),axisMax.getZ(), width, axisColor ));

		// draw screw component
		boolean positiveScrew = Math.signum(rotationAxis.getX()) == Math.signum(screwTranslation.getX());
		if( positiveScrew ) {
			// screw is in the same direction as the axis
			result.append( String.format(
					"draw ID screw"+axisID+" VECTOR {%f,%f,%f} {%f,%f,%f} WIDTH %f COLOR %s ;",
					axisMax.getX(),axisMax.getY(),axisMax.getZ(),
					screwTranslation.getX(),screwTranslation.getY(),screwTranslation.getZ(),
					width, screwColor ));
		} else {
			// screw is in the opposite direction as the axis
			result.append( String.format(
					"draw ID screw"+axisID+" VECTOR {%f,%f,%f} {%f,%f,%f} WIDTH %f COLOR %s ;",
					axisMin.getX(),axisMin.getY(),axisMin.getZ(),
					screwTranslation.getX(),screwTranslation.getY(),screwTranslation.getZ(),
					width, screwColor ));
		}

		// draw angle of rotation
		if(rotationPos != null) {
			result.append(System.getProperty("line.separator"));
			result.append(String.format("draw ID rotArc"+axisID+" ARC {%f,%f,%f} {%f,%f,%f} {0,0,0} {0,%f,%d} SCALE 500 DIAMETER %f COLOR %s;",
					axisMin.getX(),axisMin.getY(),axisMin.getZ(),
					axisMax.getX(),axisMax.getY(),axisMax.getZ(),
					Math.toDegrees(theta),
					positiveScrew ? 0 : 1 , // draw at the opposite end from the screw arrow
							width, axisColor ));
		}

		return result.toString();
	}

	/**
	 * Projects a given point onto the axis of rotation
	 * @param point
	 * @return An atom which lies on the axis, or null if the RotationAxis is purely translational
	 */
	public Atom getProjectedPoint(Atom point) {
		if(rotationPos == null) {
			// translation only
			return null;
		}

		Atom localPoint = Calc.subtract(point, rotationPos);
		double dot = Calc.scalarProduct(localPoint, rotationAxis);

		Atom localProjected = Calc.scale(rotationAxis, dot);
		Atom projected = Calc.add(localProjected, rotationPos);
		return projected;
	}

	/**
	 * Get the distance from a point to the axis of rotation
	 * @param point
	 * @return The distance to the axis, or NaN if the RotationAxis is purely translational
	 */
	public double getProjectedDistance(Atom point) {
		Atom projected = getProjectedPoint(point);
		if( projected == null) {
			// translation only
			return Double.NaN;
		}


		return Calc.getDistance(point, projected);

	}

	public void rotate(Atom[] atoms, double theta) {
		Matrix rot = getRotationMatrix(theta);
		if(rotationPos == null) {
			// Undefined rotation axis; do nothing
			return;
		}
		Atom negPos;

			negPos = Calc.invert(rotationPos);

		for(Atom a: atoms) {
			Calc.shift(a, negPos);
		}
		Calc.rotate(atoms, rot);
		for(Atom a: atoms) {
			Calc.shift(a, rotationPos);
		}
	}

	/**
	 * Calculate the rotation angle for a structure
	 * @param afpChain
	 * @return The rotation angle, in radians
	 * @throws StructureException If the alignment doesn't contain any blocks
	 * @throws NullPointerException If the alignment doesn't have a rotation matrix set
	 */
	public static double getAngle(AFPChain afpChain) throws StructureException {
		if(afpChain.getBlockNum() < 1) {
			throw new StructureException("No aligned residues");
		}
		Matrix rotation = afpChain.getBlockRotationMatrix()[0];

		if(rotation == null) {
			throw new NullPointerException("AFPChain does not contain a rotation matrix");
		}
		return getAngle(rotation);
	}
	/**
	 * Calculate the rotation angle for a given matrix
	 * @param rotation Rotation matrix
	 * @return The angle, in radians
	 */
	public static double getAngle(Matrix rotation) {
		double c = (rotation.trace()-1)/2.0; //=cos(theta)
		return Math.acos(c);
	}

	/**
	 *
	 * @return If the rotation axis is well defined, rather than purely translational
	 */
	public boolean isDefined() {
		return rotationPos != null;
	}
}
