package org.biojava.bio.structure.align.util;

import java.io.StringWriter;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.jama.Matrix;

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
 * angles. Therefor it's direction is left as null for angles less than 
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

	/**
	 * Get the component of translation perpendicular to the axis of rotation.
	 * This isn't particularly meaningful, but is calculated internally and
	 * was useful for debugging.
	 * @return
	 */
	@Deprecated
	public Atom getOtherTranslation() {
		return otherTranslation;
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
	 * Determine the location of the rotation axis based on a rotation matrix and a translation vector
	 * @param rotation
	 * @param translation
	 */
	public RotationAxis(Matrix rotation, Atom translation) {
		init(rotation, translation);
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
		double dotProduct = Calc.skalarProduct(rotationAxis, translation);

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
	 * 
	 * <p>As the rotation angle gets smaller, the axis of rotation becomes poorly
	 * defined and would need to get farther and farther away from the protein.
	 * This is not particularly useful, so we arbitrarily draw it parallel to
	 * the translation and omit the arc.
	 * @param atoms Some atoms from the protein, used for determining the bounds
	 *  	  of the axis.
	 * @return The Jmol script, suitable for calls to
	 * {@link org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol#evalString() jmol.evalString()}
	 */
	public String getJmolScript(Atom[] atoms){
		final double width=.5;// width of JMol object
		final String axisColor = "yellow"; //axis color
		final String screwColor = "orange"; //screw translation color
		
		// Project each Atom onto the rotation axis to determine limits
		double min, max;
		min = max = Calc.skalarProduct(rotationAxis,atoms[0]);
		for(int i=1;i<atoms.length;i++) {
			double prod = Calc.skalarProduct(rotationAxis,atoms[i]);
			if(prod<min) min = prod;
			if(prod>max) max = prod;
		}
		double uLen = Calc.skalarProduct(rotationAxis,rotationAxis);// Should be 1, but double check
		min/=uLen;
		max/=uLen;

		// Project the origin onto the axis. If the axis is undefined, use the center of mass
		Atom axialPt;
		if(rotationPos == null) {
			Atom center = Calc.centerOfMass(atoms);

			// Project center onto the axis
			Atom centerOnAxis = Calc.scale(rotationAxis, Calc.skalarProduct(center, rotationAxis));

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

		StringWriter result = new StringWriter();
		
		// set arrow heads to a reasonable length
		result.append("set defaultDrawArrowScale 2.0;");
		
		// draw axis of rotation
		result.append(	
				String.format("draw ID rot CYLINDER {%f,%f,%f} {%f,%f,%f} WIDTH %f COLOR %s ;",
						axisMin.getX(),axisMin.getY(),axisMin.getZ(),
						axisMax.getX(),axisMax.getY(),axisMax.getZ(), width, axisColor ));

		// draw screw component
		boolean positiveScrew = Math.signum(rotationAxis.getX()) == Math.signum(screwTranslation.getX()); 
		if( positiveScrew ) {
			// screw is in the same direction as the axis
			result.append( String.format(
					"draw ID screw VECTOR {%f,%f,%f} {%f,%f,%f} WIDTH %f COLOR %s ;",
					axisMax.getX(),axisMax.getY(),axisMax.getZ(),
					screwTranslation.getX(),screwTranslation.getY(),screwTranslation.getZ(),
					width, screwColor ));
		} else {
			// screw is in the opposite direction as the axis
			result.append( String.format(
					"draw ID screw VECTOR {%f,%f,%f} {%f,%f,%f} WIDTH %f COLOR %s ;",
					axisMin.getX(),axisMin.getY(),axisMin.getZ(),
					screwTranslation.getX(),screwTranslation.getY(),screwTranslation.getZ(),
					width, screwColor ));
		}
				
		// draw angle of rotation
		if(rotationPos != null) {
			result.append(System.getProperty("line.separator"));
			result.append(String.format("draw ID rotArc ARC {%f,%f,%f} {%f,%f,%f} {0,0,0} {0,%f,%d} SCALE 500 DIAMETER %f COLOR %s;",
					axisMin.getX(),axisMin.getY(),axisMin.getZ(),
					axisMax.getX(),axisMax.getY(),axisMax.getZ(),
					Math.toDegrees(theta),
					positiveScrew ? 0 : 1 , // draw at the opposite end from the screw arrow
					width, axisColor ));
		}

		return result.toString();
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
}
