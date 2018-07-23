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
package org.biojava.nbio.structure.symmetry.axis;

import org.biojava.nbio.structure.geometry.CalcPoint;
import org.biojava.nbio.structure.geometry.MomentsOfInertia;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.core.Rotation;
import org.biojava.nbio.structure.symmetry.core.RotationGroup;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetrySubunits;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.*;

import java.util.*;

public class RotationAxisAligner extends AxisAligner{
	
	private static final Logger logger = LoggerFactory
			.getLogger(RotationAxisAligner.class);
	
	private static final Vector3d X_AXIS = new Vector3d(1,0,0);
	private static final Vector3d Y_AXIS = new Vector3d(0,1,0);
	private static final Vector3d Z_AXIS = new Vector3d(0,0,1);

	private QuatSymmetrySubunits subunits = null;
	private RotationGroup rotationGroup = null;

	private Matrix4d transformationMatrix = new Matrix4d();
	private Matrix4d reverseTransformationMatrix = new Matrix4d();
	private Vector3d referenceVector = new Vector3d();
	private Vector3d principalRotationVector = new Vector3d();
	private Vector3d[] principalAxesOfInertia = null;
	 List<List<Integer>> alignedOrbits = null;

	private Vector3d minBoundary = new Vector3d();
	private Vector3d maxBoundary = new Vector3d();
	private double xyRadiusMax = Double.MIN_VALUE;

	boolean modified = true;

	public RotationAxisAligner(QuatSymmetryResults results) {
		this.subunits = new QuatSymmetrySubunits(results.getSubunitClusters());
		this.rotationGroup = results.getRotationGroup();

		if (subunits == null) {
			throw new IllegalArgumentException("AxisTransformation: Subunits are null");
		} else if (rotationGroup == null) {
			throw new IllegalArgumentException("AxisTransformation: RotationGroup is null");
		} else if (subunits.getSubunitCount() == 0) {
			throw new IllegalArgumentException("AxisTransformation: Subunits is empty");
		} else if (rotationGroup.getOrder() == 0) {
			throw new IllegalArgumentException("AxisTransformation: RotationGroup is empty");
		}
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getTransformation()
	 */
	@Override
	public String getSymmetry() {
		run();
		return rotationGroup.getPointGroup();
	}

	@Override
	public Matrix4d getTransformation() {
		run();
		return transformationMatrix;
	}

	@Override
	public Matrix3d getRotationMatrix() {
		run();
		Matrix3d m = new Matrix3d();
		transformationMatrix.getRotationScale(m);
		return m;
	}

	@Override
	public Matrix4d getReverseTransformation() {
		run();
		return reverseTransformationMatrix;
	}

	@Override
	public Vector3d getPrincipalRotationAxis() {
		run();
		return principalRotationVector;
	}

	@Override
	public Vector3d getRotationReferenceAxis() {
		run();
		return referenceVector;
	}

	@Override
	public Vector3d[] getPrincipalAxesOfInertia() {
		run();
		return principalAxesOfInertia;
	}

	@Override
	public Vector3d getDimension() {
		run();
		Vector3d dimension = new Vector3d();
		dimension.sub(maxBoundary, minBoundary);
		dimension.scale(0.5);
		return dimension;
	}

	/**
	 * Returns the radius for drawing the minor rotation axis in the xy-plane
	 * @return double radius in xy-plane
	 */
	@Override
	public double getRadius() {
		run();
		return xyRadiusMax;
	}

	/**
	 * Returns a transformation matrix transform polyhedra for Cn structures.
	 * The center in this matrix is the geometric center, rather then the centroid.
	 * In Cn structures those are usually not the same.
	 * @return
	 */
	@Override
	public Matrix4d getGeometicCenterTransformation() {
		run();

		Matrix4d geometricCentered = new Matrix4d(reverseTransformationMatrix);
		geometricCentered.setTranslation(new Vector3d(getGeometricCenter()));

		return geometricCentered;
	}

	/**
	 * Returns the geometric center of polyhedron. In the case of the Cn
	 * point group, the centroid and geometric center are usually not
	 * identical.
	 * @return
	 */
	@Override
	public Point3d getGeometricCenter() {
		run();

		Point3d geometricCenter = new Point3d();
		Vector3d translation = new Vector3d();
		reverseTransformationMatrix.get(translation);

		// calculate adjustment around z-axis and transform adjustment to
		//  original coordinate frame with the reverse transformation
		if (rotationGroup.getPointGroup().startsWith("C")) {
			Vector3d corr = new Vector3d(0,0, minBoundary.z+getDimension().z);
			reverseTransformationMatrix.transform(corr);
			geometricCenter.set(corr);
		}

		geometricCenter.add(translation);
		return geometricCenter;
	}

	@Override
	public Point3d getCentroid() {
		return new Point3d(subunits.getCentroid());
	}

	@Override
	public QuatSymmetrySubunits getSubunits() {
		return subunits;
	}

	public RotationGroup getRotationGroup() {
		return rotationGroup;
	}

	@Override
	public List<List<Integer>> getOrbits() {
		return alignedOrbits;
	}

	/**
	 * @return
	 */

	private void run () {
		if (modified) {
			calcPrincipalRotationVector();
			calcPrincipalAxes();
			// initial alignment with draft reference axis
			calcReferenceVector();
			calcTransformation();

			// refine ref. axis for cyclic and dihedral systems
			if ((rotationGroup.getPointGroup().startsWith("C") &&
					!rotationGroup.getPointGroup().startsWith("C2")) ||
				(rotationGroup.getPointGroup().startsWith("D") &&
							!rotationGroup.getPointGroup().startsWith("D2"))
					) {
				refineReferenceVector();
				calcTransformation();
			}
			calcReverseTransformation();
			calcBoundaries();
			if (! rotationGroup.getPointGroup().equals("Helical")) {
				calcAlignedOrbits();
			}
			modified = false;
		}
	}
	/**
	 * Returns a list of orbits (an orbit is a cyclic permutation of subunit indices that are related
	 * by a rotation around the principal rotation axis) ordered from the +z direction to the -z direction (z-depth).
	 * Within an orbit, subunit indices are ordered with a clockwise rotation around the z-axis.
	 * @return list of orbits ordered by z-depth
	 */
	private void calcAlignedOrbits() {
		Map<Double, List<Integer>> depthMap = new TreeMap<Double, List<Integer>>();
		double[] depth = getSubunitZDepth();
		alignedOrbits = calcOrbits();

		// calculate the mean depth of orbit along z-axis
		for (List<Integer> orbit: alignedOrbits) {
			// calculate the mean depth along z-axis for each orbit
			double meanDepth = 0;
			for (int subunit: orbit) {
				meanDepth += depth[subunit];
			}
			meanDepth /= orbit.size();

			if (depthMap.get(meanDepth) != null) {
				// System.out.println("Conflict in depthMap");
				meanDepth += 0.01;
			}
			depthMap.put(meanDepth, orbit);
		}

		// now fill orbits back into list ordered by depth
		alignedOrbits.clear();
		for (List<Integer> orbit: depthMap.values()) {
			// order subunit in a clockwise rotation around the z-axis
			/// starting at the 12 O-clock position (+y position)
			alignWithReferenceAxis(orbit);
			alignedOrbits.add(orbit);
		}
	}

	/**
	 * Returns an ordered list of subunit ids (orbit) in such a way that the subunit
	 * indices start at the 12 o-clock (+y axis) and proceed in a clockwise direction
	 * to the 11 o-clock position to close the "orbit".
	 *
	 * @param orbit list of subunit indices that are transformed into each other by a rotation
	 * @return list of subunit indices ordered in a clockwise direction
	 */

	private List<Integer> alignWithReferenceAxis(List<Integer> orbit) {
		int n = subunits.getSubunitCount();
		int fold = rotationGroup.getRotation(0).getFold();
		if (fold < 2) {
			return orbit;
		}
		Vector3d probe = new Vector3d();

		double dotMin1 = Double.MIN_VALUE;
		double dotMin2 = Double.MIN_VALUE;
		int index1 = 0;
		int index2 = 0;
		Vector3d Y1 = new Vector3d(0,1,0);
		Vector3d Y2 = new Vector3d(0,1,0);
		Matrix3d m = new Matrix3d();
		double angle = -2*Math.PI/fold;
		m.rotZ(0.1*angle); // add small offset, since two subunits may be equidistant to the y-axis
		m.transform(Y1);
		m.rotZ(1.1*angle);
		m.transform(Y2);
		// transform subunit centers into z-aligned position and calculate
		// width in xy direction.
		for (int i: orbit) {
			Point3d p = subunits.getCenters().get(i);
			probe.set(p);
			transformationMatrix.transform(probe);
			// find subunit that lines up with y-axis
			double dot1 = Y1.dot(probe);
			if (dot1 > dotMin1) {
				dotMin1 = dot1;
				index1 = i;
			}
			// find next subunit (rotated by one fold around z-axis - clockwise)
			double dot2 = Y2.dot(probe);
			if (dot2 > dotMin2) {
				dotMin2 = dot2;
				index2 = i;
			}
		}
//		System.out.println("Index1/2: " + index1 + " - " + index2);
//		System.out.println("Orbit0: " + orbit);
		// order subunit indices in a clockwise orientation around the z-axis
		// bring subunit into position 0
		for (int i = 0; i < n; i++) {
			if (orbit.get(0) == index1) {
				break;
			}
			Collections.rotate(orbit,1);
		}
//		System.out.println("Orbit1: " + orbit);
		// bring second subunit  onto position 1
		if (orbit.get(1) == index2) {
			return orbit;
		}
		Collections.reverse(orbit.subList(1,  orbit.size()));
		if (orbit.get(1) != index2) {
			logger.warn("alignWithReferenceAxis failed");
		}
//		System.out.println("Orbit2: " + orbit);
		return orbit;
	}



	private void calcTransformation() {
		if (rotationGroup.getPointGroup().equals("C1")) {
			calcTransformationByInertiaAxes();
		} else {
			calcTransformationBySymmetryAxes();
		}
		// make sure this value is zero. On Linux this value is 0.
		transformationMatrix.setElement(3, 3, 1.0);
	}

	private void calcReverseTransformation() {
		reverseTransformationMatrix.invert(transformationMatrix);
	}

	private void calcTransformationBySymmetryAxes() {
		Vector3d[] axisVectors = new Vector3d[2];
		axisVectors[0] = new Vector3d(principalRotationVector);
		axisVectors[1] = new Vector3d(referenceVector);

		//  y,z axis centered at the centroid of the subunits
		Vector3d[] referenceVectors = new Vector3d[2];
		referenceVectors[0] = new Vector3d(Z_AXIS);
		referenceVectors[1] = new Vector3d(Y_AXIS);

		transformationMatrix = alignAxes(axisVectors, referenceVectors);

		// combine with translation
		Matrix4d combined = new Matrix4d();
		combined.setIdentity();
		Vector3d trans = new Vector3d(subunits.getCentroid());
		trans.negate();
		combined.setTranslation(trans);
		transformationMatrix.mul(combined);

		// for cyclic geometry, set a canonical view for the Z direction
		if (rotationGroup.getPointGroup().startsWith("C")) {
			calcZDirection();
		}
	}

	private void calcTransformationByInertiaAxes() {
		Vector3d[] axisVectors = new Vector3d[2];
		axisVectors[0] = new Vector3d(principalAxesOfInertia[0]);
		axisVectors[1] = new Vector3d(principalAxesOfInertia[1]);

		Vector3d[] referenceVectors = new Vector3d[2];
		referenceVectors[0] = new Vector3d(Y_AXIS);
		referenceVectors[1] = new Vector3d(X_AXIS);

		// align inertia axes with y-x plane
		transformationMatrix = alignAxes(axisVectors, referenceVectors);

		// combine with translation
		Matrix4d translation = new Matrix4d();
		translation.setIdentity();
		Vector3d trans = new Vector3d(subunits.getCentroid());
		trans.negate();
		translation.setTranslation(trans);
		transformationMatrix.mul(translation);
	}

	/**
	 * Returns a transformation matrix that rotates refPoints to match
	 * coordPoints
	 * @param refPoints the points to be aligned
	 * @param referenceVectors
	 * @return
	 */
	private Matrix4d alignAxes(Vector3d[] axisVectors, Vector3d[] referenceVectors) {
		Matrix4d m1 = new Matrix4d();
		AxisAngle4d a = new AxisAngle4d();
		Vector3d axis = new Vector3d();

		// calculate rotation matrix to rotate refPoints[0] into coordPoints[0]
		Vector3d v1 = new Vector3d(axisVectors[0]);
		Vector3d v2 = new Vector3d(referenceVectors[0]);
		double dot = v1.dot(v2);
		if (Math.abs(dot) < 0.999) {
			axis.cross(v1,v2);
			axis.normalize();
			a.set(axis, v1.angle(v2));
			m1.set(a);
			// make sure matrix element m33 is 1.0. It's 0 on Linux.
			m1.setElement(3,  3, 1.0);
		} else if (dot > 0) {
			// parallel axis, nothing to do -> identity matrix
			m1.setIdentity();
		} else if (dot < 0) {
			// anti-parallel axis, flip around x-axis
			m1.set(flipX());
		}

		// apply transformation matrix to all refPoints
		m1.transform(axisVectors[0]);
		m1.transform(axisVectors[1]);

		// calculate rotation matrix to rotate refPoints[1] into coordPoints[1]
		v1 = new Vector3d(axisVectors[1]);
		v2 = new Vector3d(referenceVectors[1]);
		Matrix4d m2 = new Matrix4d();
		dot = v1.dot(v2);
		if (Math.abs(dot) < 0.999) {
			axis.cross(v1,v2);
			axis.normalize();
			a.set(axis, v1.angle(v2));
			m2.set(a);
			// make sure matrix element m33 is 1.0. It's 0 on Linux.
			m2.setElement(3,  3, 1.0);
		} else if (dot > 0) {
			// parallel axis, nothing to do -> identity matrix
			m2.setIdentity();
		} else if (dot < 0) {
			// anti-parallel axis, flip around z-axis
			m2.set(flipZ());
		}

		// apply transformation matrix to all refPoints
		m2.transform(axisVectors[0]);
		m2.transform(axisVectors[1]);

		// combine the two rotation matrices
		m2.mul(m1);

		// the RMSD should be close to zero
		Point3d[] axes = new Point3d[2];
		axes[0] = new Point3d(axisVectors[0]);
		axes[1] = new Point3d(axisVectors[1]);
		Point3d[] ref = new Point3d[2];
		ref[0] = new Point3d(referenceVectors[0]);
		ref[1] = new Point3d(referenceVectors[1]);
		if (CalcPoint.rmsd(axes, ref) > 0.1) {
			logger.warn("AxisTransformation: axes alignment is off. RMSD: " 
					+ CalcPoint.rmsd(axes, ref));
		}

		return m2;
	}

	private void calcPrincipalAxes() {
		MomentsOfInertia moi = new MomentsOfInertia();

		for (Point3d[] list: subunits.getTraces()) {
			for (Point3d p: list) {
				moi.addPoint(p, 1.0);
			}
		}
		principalAxesOfInertia = moi.getPrincipalAxes();
	}

	/**
	 * Calculates the min and max boundaries of the structure after it has been
	 * transformed into its canonical orientation.
	 */
	private void calcBoundaries() {
		minBoundary.x = Double.MAX_VALUE;
		maxBoundary.x = Double.MIN_VALUE;
		minBoundary.y = Double.MAX_VALUE;
		maxBoundary.x = Double.MIN_VALUE;
		minBoundary.z = Double.MAX_VALUE;
		maxBoundary.z = Double.MIN_VALUE;

		Point3d probe = new Point3d();

		for (Point3d[] list: subunits.getTraces()) {
			for (Point3d p: list) {
				probe.set(p);
				transformationMatrix.transform(probe);

				minBoundary.x = Math.min(minBoundary.x, probe.x);
				maxBoundary.x = Math.max(maxBoundary.x, probe.x);
				minBoundary.y = Math.min(minBoundary.y, probe.y);
				maxBoundary.y = Math.max(maxBoundary.y, probe.y);
				minBoundary.z = Math.min(minBoundary.z, probe.z);
				maxBoundary.z = Math.max(maxBoundary.z, probe.z);
				xyRadiusMax = Math.max(xyRadiusMax, Math.sqrt(probe.x*probe.x + probe.y * probe.y));
			}
		}
	}

	/*
	 * Modifies the rotation part of the transformation axis for
	 * a Cn symmetric complex, so that the narrower end faces the
	 * viewer, and the wider end faces away from the viewer. Example: 3LSV
	 */
	private void calcZDirection() {
		calcBoundaries();

		// if the longer part of the structure faces towards the back (-z direction),
		// rotate around y-axis so the longer part faces the viewer (+z direction)
		if (Math.abs(minBoundary.z) > Math.abs(maxBoundary.z)) {
			Matrix4d rot = flipY();
			rot.mul(transformationMatrix);
			transformationMatrix.set(rot);
		}
	}

	/**
	 *
	 */
	private List<List<Integer>> getOrbitsByXYWidth() {
		Map<Double, List<Integer>> widthMap = new TreeMap<Double, List<Integer>>();
		double[] width = getSubunitXYWidth();
		List<List<Integer>> orbits = calcOrbits();

		// calculate the mean width of orbits in XY-plane
		for (List<Integer> orbit: orbits) {
			double meanWidth = 0;
			for (int subunit: orbit) {
				meanWidth += width[subunit];
			}
			meanWidth /= orbit.size();

			if (widthMap.get(meanWidth) != null) {
				meanWidth += 0.01;
			}
			widthMap.put(meanWidth, orbit);
		}

		// now fill orbits back into list ordered by width
		orbits.clear();
		for (List<Integer> orbit: widthMap.values()) {
			orbits.add(orbit);
		}
		return orbits;
	}

	private double[] getSubunitXYWidth() {
		int n = subunits.getSubunitCount();
		double[] width = new double[n];
		Point3d probe = new Point3d();

		// transform subunit centers into z-aligned position and calculate
		// width in xy direction.
		for (int i = 0; i < n; i++) {
			width[i] = Double.MIN_VALUE;
			for (Point3d p: subunits.getTraces().get(i)) {
				probe.set(p);
				transformationMatrix.transform(probe);
				width[i] = Math.max(width[i], Math.sqrt(probe.x*probe.x + probe.y*probe.y));
			}
		}
		return width;
	}

	private double[] getSubunitZDepth() {
		int n = subunits.getSubunitCount();
		double[] depth = new double[n];
		Point3d probe = new Point3d();

		// transform subunit centers into z-aligned position and calculate
		// z-coordinates (depth) along the z-axis.
		for (int i = 0; i < n; i++) {
			Point3d p= subunits.getCenters().get(i);
			probe.set(p);
			transformationMatrix.transform(probe);
			depth[i] = probe.z;
		}
		return depth;
	}

	/**
	 * Returns a list of list of subunit ids that form an "orbit", i.e. they
	 * are transformed into each other during a rotation around the principal symmetry axis (z-axis)
	 * @return
	 */
	private List<List<Integer>> calcOrbits() {
		int n = subunits.getSubunitCount();
		int fold = rotationGroup.getRotation(0).getFold();

		List<List<Integer>> orbits = new ArrayList<List<Integer>>();
		boolean[] used = new boolean[n];
		Arrays.fill(used, false);

		for (int i = 0; i < n; i++) {
			if (! used[i]) {
				// determine the equivalent subunits
				List<Integer> orbit = new ArrayList<Integer>(fold);
				for (int j = 0; j < fold; j++) {
					List<Integer> permutation = rotationGroup.getRotation(j).getPermutation();
					orbit.add(permutation.get(i));
					used[permutation.get(i)] = true;
				}
				orbits.add(deconvolute(orbit));
			}
		}
		return orbits;
	}

	private List<Integer> deconvolute(List<Integer> orbit) {
		if (rotationGroup.getOrder() < 2) {
			return orbit;
		}
		List<Integer> p0 = rotationGroup.getRotation(0).getPermutation();
		List<Integer> p1 = rotationGroup.getRotation(1).getPermutation();
//		System.out.println("deconvolute");
//		System.out.println("Permutation0: " + p0);
//		System.out.println("Permutation1: " + p1);

		List<Integer> inRotationOrder = new ArrayList<Integer>(orbit.size());
		inRotationOrder.add(orbit.get(0));
		for (int i = 1; i < orbit.size(); i++) {
			int current = inRotationOrder.get(i-1);
			int index = p0.indexOf(current);
			int next = p1.get(index);
			if (!orbit.contains(next)) {
				logger.warn("deconvolute: inconsistency in permuation. Returning original order");
				return orbit;
			}
			inRotationOrder.add(next);
		}
//		System.out.println("In order: " + inRotationOrder);
		return inRotationOrder;
	}

	/**
	 * Returns a vector along the principal rotation axis for the
	 * alignment of structures along the z-axis
	 * @return principal rotation vector
	 */
	private void calcPrincipalRotationVector() {
		Rotation rotation = rotationGroup.getRotation(0); // the rotation around the principal axis is the first rotation
		AxisAngle4d axisAngle = rotation.getAxisAngle();
		principalRotationVector = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
	}

	/**
	 * Returns a vector perpendicular to the principal rotation vector
	 * for the alignment of structures in the xy-plane
	 * @return reference vector
	 */
	private void calcReferenceVector() {
		referenceVector = null;
		if (rotationGroup.getPointGroup().startsWith("C")) {
			referenceVector = getReferenceAxisCylic();
		} else if (rotationGroup.getPointGroup().startsWith("D")) {
			referenceVector = getReferenceAxisDihedral();
		} else if (rotationGroup.getPointGroup().equals("T")) {
			referenceVector = getReferenceAxisTetrahedral();
		} else if (rotationGroup.getPointGroup().equals("O")) {
			referenceVector = getReferenceAxisOctahedral();
		} else if (rotationGroup.getPointGroup().equals("I")) {
			referenceVector = getReferenceAxisIcosahedral();
		} else if (rotationGroup.getPointGroup().equals("Helical")) {
			// TODO what should the reference vector be??
			referenceVector = getReferenceAxisCylic();
		}

		if (referenceVector == null) {
			logger.warn("no reference vector found. Using y-axis.");
			referenceVector = new Vector3d(Y_AXIS);
		}
		// make sure reference vector is perpendicular principal roation vector
		referenceVector = orthogonalize(principalRotationVector, referenceVector);
	}

	/**
	 * Returns a normalized vector that represents a minor rotation axis, except
	 * for Cn, this represents an axis orthogonal to the principal axis.
	 * @return minor rotation axis
	 */
	private void refineReferenceVector() {
		referenceVector = new Vector3d(Y_AXIS);
		if (rotationGroup.getPointGroup().startsWith("C")) {
			referenceVector = getReferenceAxisCylicWithSubunitAlignment();
		} else if (rotationGroup.getPointGroup().startsWith("D")) {
			referenceVector = getReferenceAxisDihedralWithSubunitAlignment();
		}

		referenceVector = orthogonalize(principalRotationVector, referenceVector);
	}

	private Vector3d orthogonalize(Vector3d vector1, Vector3d vector2) {
		double dot = vector1.dot(vector2);
		Vector3d ref = new Vector3d(vector2);
//		System.out.println("p.r: " + dot);
//		System.out.println("Orig refVector: " + referenceVector);
		if (dot < 0) {
			vector2.negate();
		}
		vector2.cross(vector1, vector2);
//		System.out.println("Intermed. refVector: " + vector2);
		vector2.normalize();
//		referenceVector.cross(referenceVector, principalRotationVector);
		vector2.cross(vector1, vector2);
		vector2.normalize();
		if (ref.dot(vector2) < 0) {
			vector2.negate();
		}
//		System.out.println("Mod. refVector: " + vector2);
		return vector2;
	}
	/**
	 * Returns the default reference vector for the alignment of Cn structures
	 * @return
	 */
	private Vector3d getReferenceAxisCylic() {
		if (rotationGroup.getPointGroup().equals("C2")) {
			Vector3d vr = new Vector3d(subunits.getOriginalCenters().get(0));
			vr.sub(subunits.getCentroid());
			vr.normalize();
			return vr;
		}

		// get principal axis vector that is perpendicular to the principal
		// rotation vector
		Vector3d vmin = null;
		double dotMin = 1.0;
		for (Vector3d v: principalAxesOfInertia) {
			if (Math.abs(principalRotationVector.dot(v)) < dotMin) {
				dotMin = Math.abs(principalRotationVector.dot(v));
				vmin = new Vector3d(v);
			}
		}
		if (principalRotationVector.dot(vmin) < 0) {
			vmin.negate();
		}

		return vmin;
	}


	/**
	 * Returns a reference vector for the alignment of Cn structures.
	 * @return reference vector
	 */
	private Vector3d getReferenceAxisCylicWithSubunitAlignment() {
		if (rotationGroup.getPointGroup().equals("C2")) {
			return referenceVector;
		}

		// find subunit that extends the most in the xy-plane
		List<List<Integer>> orbits = getOrbitsByXYWidth();
		// get the last orbit which is the widest
		List<Integer> widestOrbit = orbits.get(orbits.size()-1);
		List<Point3d> centers = subunits.getCenters();
		int subunit = widestOrbit.get(0);

		// calculate reference vector
		Vector3d refAxis = new Vector3d();
		refAxis.sub(centers.get(subunit), subunits.getCentroid());
		refAxis.normalize();
		return refAxis;
	}

	/**
	 *
	 */
	private Vector3d getReferenceAxisDihedralWithSubunitAlignment() {
		int maxFold = rotationGroup.getRotation(0).getFold();

		double minAngle = Double.MAX_VALUE;
		Vector3d refVec = null;

		Vector3d ref = getSubunitReferenceVector();

		for (int i = 0; i < rotationGroup.getOrder(); i++) {
			if (rotationGroup.getRotation(i).getDirection() == 1 &&
					(rotationGroup.getRotation(i).getFold() < maxFold) ||
					rotationGroup.getPointGroup().equals("D2")) {

				AxisAngle4d axisAngle = rotationGroup.getRotation(i).getAxisAngle();
				Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				v.normalize();

//				System.out.println("Ref axis angle(+): " + Math.toDegrees(v.angle(ref)));
				double angle =  v.angle(ref);
				if (angle < minAngle) {
					minAngle = angle;
					refVec = v;
				}
				Vector3d vn = new Vector3d(v);
				vn.negate();
//				System.out.println("Ref axis angle(-): " + Math.toDegrees(vn.angle(ref)));
				angle =  vn.angle(ref);
				if (angle < minAngle) {
					minAngle = angle;
					refVec = vn;
				}
			}
		}
		refVec.normalize();
		return refVec;
	}

	/**
	 *
	 */
	private Vector3d getReferenceAxisDihedral() {
		int maxFold = rotationGroup.getRotation(0).getFold();
		// one exception: D2
		if (maxFold == 2) {
			maxFold = 3;
		}
		// TODO how about D2, where minor axis = 2 = principal axis??
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
			if (rotationGroup.getRotation(i).getDirection() == 1 && rotationGroup.getRotation(i).getFold() < maxFold) {
				AxisAngle4d axisAngle = rotationGroup.getRotation(i).getAxisAngle();
				Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				v.normalize();
				return v;
			}
		}
		return null;
	}

	private Vector3d getReferenceAxisTetrahedral() {
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
				AxisAngle4d axisAngle = rotationGroup.getRotation(i).getAxisAngle();
				Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				double d = v.dot(principalRotationVector);
				if (rotationGroup.getRotation(i).getFold() == 3) {
					// the dot product 0 is between to adjacent 3-fold axes
					if (d > 0.3 && d < 0.9) {
						return v;
					}
				}
		}
		return null;
	}

	private Vector3d getReferenceAxisOctahedral() {
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
				AxisAngle4d axisAngle = rotationGroup.getRotation(i).getAxisAngle();
				Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				double d = v.dot(principalRotationVector);
				if (rotationGroup.getRotation(i).getFold() == 4) {
					// the dot product 0 is between to adjacent 4-fold axes
					if (d > -0.1 && d < 0.1 ) {
						return v;
					}
				}
		}
		return null;
	}

	private Vector3d getReferenceAxisIcosahedral() {
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
				AxisAngle4d axisAngle = rotationGroup.getRotation(i).getAxisAngle();
				Vector3d v = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				double d = v.dot(principalRotationVector);
				if (rotationGroup.getRotation(i).getFold() == 5) {
					// the dot product of 0.447.. is between to adjacent 5-fold axes
//					if (d > 0.447 && d < 0.448) {
					if (d > 0.4 && d < 0.5) {
						return v;
					}
				}
		}
		return null;
	}

	private Vector3d getSubunitReferenceVector() {
			int n = subunits.getSubunitCount();
			Point3d probe = new Point3d();

			// transform subunit centers into z-aligned position and calculate
			// width in xy direction.
			double maxWidthSq = 0;
			Point3d ref = null;
			for (int i = 0; i < n; i++) {
				for (Point3d p: subunits.getTraces().get(i)) {
					probe.set(p);
					transformationMatrix.transform(probe);
					double widthSq = probe.x*probe.x + probe.y*probe.y;
					if (widthSq > maxWidthSq) {
						maxWidthSq = widthSq;
						ref = p;
					}
				}
			}
	//		System.out.println("width: " + maxWidthSq);
			Vector3d refVector = new Vector3d();
			refVector.sub(ref, subunits.getCentroid());
			refVector.normalize();
			return refVector;
		}

	private static Matrix4d flipX() {
		Matrix4d rot = new Matrix4d();
		rot.m00 = 1;
		rot.m11 = -1;
		rot.m22 = -1;
		rot.m33 = 1;
		return rot;
	}

	private static Matrix4d flipY() {
		Matrix4d rot = new Matrix4d();
		rot.m00 = -1;
		rot.m11 = 1;
		rot.m22 = -1;
		rot.m33 = 1;
		return rot;
	}

	private static Matrix4d flipZ() {
		Matrix4d rot = new Matrix4d();
		rot.m00 = -1;
		rot.m11 = -1;
		rot.m22 = 1;
		rot.m33 = 1;
		return rot;
	}

}
