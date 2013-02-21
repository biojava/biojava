package org.biojava.bio.structure.gui;

import java.io.IOException;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;

/**
 * Calculates the rotation axis for an alignment
 * @author Spencer Bliven
 *
 */
public final class RotationAxis {
	static final double MIN_ANGLE = Math.toRadians(5.);

	private final double theta;
	private final Atom rotationAxis; // axis of rotation
	private final Atom rotationPos; // a point on the axis of rotation
	private final Atom screwTranslation; //translation parallel to the axis of rotation
	private final Atom otherTranslation; // translation perpendicular to the axis of rotation
	
	/**
	 * The rotation angle
	 * @return theta
	 */
	public double getTheta() {
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
	 * Get the component of translation perpendicular to the axis of rotation
	 * @return
	 */
	public Atom getOtherTranslation() {
		return otherTranslation;
	}

	/**
	 * Determine the location of the rotation axis based on a rotation matrix and a translation vector
	 * @param rotation
	 * @param translation
	 */
	public RotationAxis(Matrix rotation, Atom translation) {
		if(rotation.getColumnDimension() != 3 || rotation.getRowDimension() != 3) {
			throw new IllegalArgumentException("Expected 3x3 rotation matrix");
		}
		
		
		// Calculate angle
		double c = (rotation.trace()-1)/2.0; //=cos(theta)
		this.theta = Math.acos(c);
		
		if(theta < MIN_ANGLE) {
			//TODO
			throw new UnsupportedOperationException("Small angles not implemented");
		}
		
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
		if( Math.abs(d0) < Math.abs(d1) ) {
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
		} else { //d0
			if(d0>=0) { //u[0] positive
				if( s01 < 0 ) rotAx[0] = -rotAx[0];
				if( s02 < 0 ) rotAx[2] = -rotAx[2];
			} else { //u[0] negative
				rotAx[0] = -rotAx[0];
				if( s01 >= 0 ) rotAx[0] = -rotAx[0];
				if( s02 >= 0 ) rotAx[2] = -rotAx[2];
			}
		}

		rotationAxis = new AtomImpl();
		rotationAxis.setCoords(rotAx);
		
		// Calculate screw = (rotationAxis dot translation)*u
		double dotProduct = Calc.skalarProduct(rotationAxis, translation);
		
		screwTranslation = Calc.scale(rotationAxis, dotProduct);
		otherTranslation = Calc.subtract(translation, screwTranslation);

		Atom hypot = Calc.vectorProduct(rotationAxis, otherTranslation);
		Calc.scaleEquals(hypot,.5/Math.tan(theta/2.0));
		
		// Calculate rotation axis position
		rotationPos = Calc.scaleAdd(.5,otherTranslation, hypot);
		
	}
	
	public void displayRotationAxis(StructureAlignmentJmol jmolPanel, Atom[] atoms) {
		
		// Project each Atom onto the rotation axis to determine limits
		double min, max, mean;
		min = max = mean = Calc.skalarProduct(rotationAxis,atoms[0]);
		for(int i=1;i<atoms.length;i++) {
			double prod = Calc.skalarProduct(rotationAxis,atoms[i]);
			mean += prod;
			if(prod<min) min = prod;
			if(prod>max) max = prod;
		}
		double uLen = Calc.skalarProduct(rotationAxis,rotationAxis);
		min/=uLen;
		max/=uLen;
		mean/=atoms.length;
		
		Atom axisMin = (Atom) rotationPos.clone();
		Calc.scaleAdd(min, rotationAxis, axisMin);
		Atom axisMax = (Atom) rotationPos.clone();
		Calc.scaleAdd(max, rotationAxis, axisMax);
		Atom axisMean = (Atom) rotationPos.clone();
		Calc.scaleAdd(mean,rotationAxis, axisMean);
		
		// draw axis
		double width=.5;
		jmolPanel.evalString(String.format("draw ID rot ARROW {%f,%f,%f} {%f,%f,%f} WIDTH %f COLOR cyan ;",
				axisMin.getX(),axisMin.getY(),axisMin.getZ(),
				axisMax.getX(),axisMax.getY(),axisMax.getZ(), width ));

		jmolPanel.evalString(String.format("draw ID rotArc ARROW ARC {%f,%f,%f} {%f,%f,%f} {0,0,0} {0,%f,1} SCALE 500 DIAMETER %f COLOR cyan;",
				axisMin.getX(),axisMin.getY(),axisMin.getZ(),
				axisMax.getX(),axisMax.getY(),axisMax.getZ(),
				Math.toDegrees(theta), width ));
	}
	
	public static void main(String[] args) {
		
		// Compare two chains of a dimer to force CE to give a symmetric alignment.
		String name1 = "1AVD.A";
		String name2 = "1AVD.B";
		String display = "1AVD";
		
//		name1 = "4HHB.A:,B:";
//		name2 = "4HHB.C:,D:";
//		display = "4HHB";
		
		AtomCache cache = new AtomCache();
		try {
			Atom[] ca1 = cache.getAtoms(name1);
			Atom[] ca2 = cache.getAtoms(name2);
			Atom[] caD = cache.getAtoms(display);
			
			StructureAlignment ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			
			AFPChain afpChain = ce.align(ca1, ca2);
			Matrix mat = afpChain.getBlockRotationMatrix()[0];
			Atom shift = afpChain.getBlockShiftVector()[0];
			
			System.out.println("Shift:");
			System.out.println(shift);
			System.out.println("Matrix:");
			System.out.println(mat);
			
			RotationAxis axis = new RotationAxis(mat,shift);

			double theta = Math.toDegrees(axis.getTheta());
			System.out.format("Angle: %f degrees%n",theta);
			
			//StructureAlignmentJmol jmolPanel = StructureAlignmentDisplay.display(afpChain, caD, caD);
			StructureAlignmentJmol jmolPanel = new StructureAlignmentJmol();
			jmolPanel.setAtoms(caD);
			
			jmolPanel.evalString("select * ; color chain;");
			jmolPanel.evalString("select nucleic; cartoon on;");
			jmolPanel.evalString("select *; spacefill off; wireframe off; cartoon on;  ");
			
			// draw coordinate axes
			jmolPanel.evalString("draw ID x VECTOR {0,0,0} {5,0,0} WIDTH 0.5 COLOR red \">x\";");
			jmolPanel.evalString("draw ID y VECTOR {0,0,0} {0,5,0} WIDTH 0.5 COLOR green \">y\";");
			jmolPanel.evalString("draw ID z VECTOR {0,0,0} {0,0,5} WIDTH 0.5 COLOR blue \">z\";");
			
			// draw axis
			axis.displayRotationAxis(jmolPanel, caD);
	
			/*
			// draw intermediate vectors for debugging
			jmolPanel.evalString(String.format("draw ID s VECTOR {0,0,0} {%f,%f,%f} WIDTH %f COLOR orange \">s\";",
					s.getX(),s.getY(),s.getZ(), width ));

			Atom perp = axis.getOtherTranslation();
			Atom screw = axis.getScrewTranslation();

			double uScale = 10;
			jmolPanel.evalString(String.format("draw ID u VECTOR {0,0,0} {%f,%f,%f} WIDTH %f COLOR orange \">u\";",
					uScale*u.getX(),uScale*u.getY(),uScale*u.getZ(), width ));
			
			jmolPanel.evalString(String.format("draw ID perp VECTOR {0,0,0} {%f,%f,%f} WIDTH %f COLOR yellow \">tPerp\";",
					perp.getX(),perp.getY(),perp.getZ(), width));
			jmolPanel.evalString(String.format("draw ID screw VECTOR {0,0,0} {%f,%f,%f} WIDTH %f COLOR yellow \">screw\";",
					screw.getX(),screw.getY(),screw.getZ(), width));

			jmolPanel.evalString(String.format("draw ID t VECTOR {0,0,0} {%f,%f,%f} WIDTH %f COLOR yellow \">t\";",
					shift.getX(),shift.getY(),shift.getZ(), width));
			*/
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (StructureException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
