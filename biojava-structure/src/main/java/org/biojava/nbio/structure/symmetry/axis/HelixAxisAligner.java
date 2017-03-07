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
import org.biojava.nbio.structure.symmetry.core.Helix;
import org.biojava.nbio.structure.symmetry.core.HelixLayers;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetrySubunits;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.*;

import java.util.*;

public class HelixAxisAligner extends AxisAligner {
	
	private static final Logger logger = LoggerFactory
			.getLogger(HelixAxisAligner.class);

	
	private static final Vector3d Y_AXIS = new Vector3d(0,1,0);
	private static final Vector3d Z_AXIS = new Vector3d(0,0,1);

	private QuatSymmetrySubunits subunits = null;
	private HelixLayers helixLayers = null;

	private Matrix4d transformationMatrix = new Matrix4d();
	private Matrix4d reverseTransformationMatrix = new Matrix4d();
	private Vector3d referenceVector = new Vector3d();
	private Vector3d principalRotationVector = new Vector3d();
	private Vector3d[] principalAxesOfInertia = null;
	private List<List<Integer>> alignedOrbits = null;

	private Vector3d minBoundary = new Vector3d();
	private Vector3d maxBoundary = new Vector3d();
	private double xzRadiusMax = Double.MIN_VALUE;

	boolean modified = true;

	public HelixAxisAligner(QuatSymmetryResults results) {
		this.subunits = new QuatSymmetrySubunits(results.getSubunitClusters());
		this.helixLayers = results.getHelixLayers();
		if (subunits == null) {
			throw new IllegalArgumentException("HelixAxisTransformation: Subunits are null");
		} else if (helixLayers == null) {
			throw new IllegalArgumentException("HelixAxisTransformation: HelixLayers is null");
		} else if (subunits.getSubunitCount() == 0) {
			throw new IllegalArgumentException("HelixAxisTransformation: Subunits is empty");
		} else if (helixLayers.size() == 0) {
			throw new IllegalArgumentException("HelixAxisTransformation: HelixLayers is empty");
		}
	}


	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getTransformation()
	 */
	@Override
	public String getSymmetry() {
		run();
		return "H";
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getTransformation()
	 */
	@Override
	public Matrix4d getTransformation() {
		run();
		return transformationMatrix;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getRotationMatrix()
	 */
	@Override
	public Matrix3d getRotationMatrix() {
		run();
		Matrix3d m = new Matrix3d();
		transformationMatrix.getRotationScale(m);
		return m;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getReverseTransformation()
	 */
	@Override
	public Matrix4d getReverseTransformation() {
		run();
		return reverseTransformationMatrix;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getPrincipalRotationAxis()
	 */
	@Override
	public Vector3d getPrincipalRotationAxis() {
		run();
		return principalRotationVector;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getRotationReferenceAxis()
	 */
	@Override
	public Vector3d getRotationReferenceAxis() {
		run();
		return referenceVector;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getPrincipalAxesOfInertia()
	 */
	@Override
	public Vector3d[] getPrincipalAxesOfInertia() {
		run();
		return principalAxesOfInertia;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getDimension()
	 */
	@Override
	public Vector3d getDimension() {
		run();
		Vector3d dimension = new Vector3d();
		dimension.sub(maxBoundary, minBoundary);
		dimension.scale(0.5);
		return dimension;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getXYRadius()
	 */
	@Override
	public double getRadius() {
		run();
		return xzRadiusMax;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getGeometicCenterTransformation()
	 */
	@Override
	public Matrix4d getGeometicCenterTransformation() {
		run();

		Matrix4d geometricCentered = new Matrix4d(reverseTransformationMatrix);
		geometricCentered.setTranslation(new Vector3d(getGeometricCenter()));

		return geometricCentered;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getGeometricCenter()
	 */
	@Override
	public Point3d getGeometricCenter() {
		run();

		Point3d geometricCenter = new Point3d();
		Vector3d translation = new Vector3d();
//		reverseTransformationMatrix.get(translation);

		// TODO does this apply to the helic case?
		// calculate adjustment around z-axis and transform adjustment to
		//  original coordinate frame with the reverse transformation

//		Vector3d corr = new Vector3d(0,minBoundary.y+getDimension().y, 0);
//		reverseTransformationMatrix.transform(corr);
//		geometricCenter.set(corr);

		reverseTransformationMatrix.transform(translation);
		geometricCenter.add(translation);
		return geometricCenter;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getCentroid()
	 */
	@Override
	public Point3d getCentroid() {
		return new Point3d(subunits.getCentroid());
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getSubunits()
	 */
	@Override
	public QuatSymmetrySubunits getSubunits() {
		return subunits;
	}

	public HelixLayers getHelixLayers() {
		run();
		return helixLayers;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.core.AxisAligner#getOrbits()
	 */
	@Override
	public List<List<Integer>> getOrbits() {
		run();
		return alignedOrbits;
	}

	/**
	 * @return
	 */

	private void run () {
		if (modified) {
			calcPrincipalRotationVector();
			calcPrincipalAxes();
			calcReferenceVector();
			calcTransformation();
			// orient helix along Y axis by rotating 90 degrees around X-axis
			transformationMatrix = reorientHelix(0);

			calcReverseTransformation();
			calcBoundaries();
			calcAlignedOrbits();
			calcCenterOfRotation();

			// orient helix along Y axis by rotating 90 degrees around X-axis
//			transformationMatrix = reorientHelix(0);
//			calcReverseTransformation();

			modified = false;
		}
	}

	public Point3d calcCenterOfRotation() {
		List<Integer> line = getLongestLayerLine();

		// can't determine center of rotation if there are only 2 points
		// TODO does this ever happen??
		if (line.size() < 3) {
			return subunits.getCentroid();
		}

		Point3d centerOfRotation = new Point3d();
		List<Point3d> centers = subunits.getOriginalCenters();

		// calculate helix mid points for each set of 3 adjacent subunits
		for (int i = 0; i < line.size()-2; i++) {
			Point3d p1 = new Point3d(centers.get(line.get(i)));
			Point3d p2 = new Point3d(centers.get(line.get(i+1)));
			Point3d p3 = new Point3d(centers.get(line.get(i+2)));
			transformationMatrix.transform(p1);
			transformationMatrix.transform(p2);
			transformationMatrix.transform(p3);
			centerOfRotation.add(getMidPoint(p1, p2, p3));
		}

		// average over all midpoints to find best center of rotation
		centerOfRotation.scale(1/(line.size()-2));
		// since helix is aligned along the y-axis, with an origin at y = 0, place the center of rotation there
		centerOfRotation.y = 0;
		// transform center of rotation to the original coordinate frame
		reverseTransformationMatrix.transform(centerOfRotation);
//		System.out.println("center of rotation: " + centerOfRotation);
		return centerOfRotation;
	}

	private List<Integer> getLongestLayerLine() {
		int len = 0;
		int index = 0;

		Helix helix = helixLayers.getByLargestContacts();
		List<List<Integer>> layerLines = helix.getLayerLines();
		for (int i = 0; i < layerLines.size(); i++) {
			if (layerLines.get(i).size() > len) {
				len = layerLines.get(i).size();
				index = i;
			}
		}
		return layerLines.get(index);
	}

	/**
	 * Return a midpoint of a helix, calculated from three positions
	 * of three adjacent subunit centers.
	 * @param p1 center of first subunit
	 * @param p2 center of second subunit
	 * @param p3 center of third subunit
	 * @return midpoint of helix
	 */
	private Point3d getMidPoint(Point3d p1, Point3d p2, Point3d p3) {
		Vector3d v1 = new Vector3d();
		v1.sub(p1, p2);
		Vector3d v2 = new Vector3d();
		v2.sub(p3, p2);
		Vector3d v3 = new Vector3d();
		v3.add(v1);
		v3.add(v2);
		v3.normalize();

		// calculat the total distance between to subunits
		double dTotal = v1.length();
		// calculate the rise along the y-axis. The helix axis is aligned with y-axis,
		// therfore, the rise between subunits is the y-distance
		double rise = p2.y - p1.y;
		// use phythagorean theoremm to calculate chord length between two subunit centers
		double chord = Math.sqrt(dTotal*dTotal - rise*rise);
//		System.out.println("Chord d: " + dTotal + " rise: " + rise + "chord: " + chord);
		double angle = helixLayers.getByLargestContacts().getAxisAngle().angle;

		// using the axis angle and the chord length, we can calculate the radius of the helix
		// http://en.wikipedia.org/wiki/Chord_%28geometry%29
		double radius = chord/Math.sin(angle/2)/2; // can this go to zero?
//		System.out.println("Radius: " + radius);

		// project the radius onto the vector that points toward the helix axis
		v3.scale(radius);
		v3.add(p2);
//		System.out.println("Angle: " + Math.toDegrees(helixLayers.getByLowestAngle().getAxisAngle().angle));
		Point3d cor = new Point3d(v3);
		return cor;
	}

	private Matrix4d reorientHelix(int index) {
		Matrix4d matrix = new Matrix4d();
		matrix.setIdentity();
		matrix.setRotation(new AxisAngle4d(1,0,0,Math.PI/2*(index+1)));
		matrix.mul(transformationMatrix);
	return matrix;
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
				meanDepth += 0.01;
			}
//			System.out.println("calcAlignedOrbits: " + meanDepth + " orbit: " + orbit);
			depthMap.put(meanDepth, orbit);
		}

		// now fill orbits back into list ordered by depth
		alignedOrbits.clear();
		for (List<Integer> orbit: depthMap.values()) {
			// order subunit in a clockwise rotation around the z-axis
			/// starting at the 12 O-clock position (+y position)
			// TODO how should this be aligned??
	//		alignWithReferenceAxis(orbit);
			alignedOrbits.add(orbit);
		}
	}


	private void calcTransformation() {
		calcTransformationBySymmetryAxes();
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

		// for helical geometry, set a canonical view for the Z direction
		calcZDirection();
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
		xzRadiusMax = Double.MIN_VALUE;

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
				xzRadiusMax = Math.max(xzRadiusMax, Math.sqrt(probe.x*probe.x + probe.z * probe.z));
			}
		}
//		System.out.println("MinBoundary: " + minBoundary);
//		System.out.println("MaxBoundary: " + maxBoundary);
//		System.out.println("zxRadius: " + xzRadiusMax);
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

		List<List<Integer>> orbits = new ArrayList<List<Integer>>();
		for (int i = 0; i < n; i++) {
			orbits.add(Collections.singletonList(i));
		}

		return orbits;
	}

	/**
	 * Returns a vector along the principal rotation axis for the
	 * alignment of structures along the z-axis
	 * @return principal rotation vector
	 */
	private void calcPrincipalRotationVector() {
//		AxisAngle4d axisAngle = helixLayers.getByLowestAngle().getAxisAngle();
		AxisAngle4d axisAngle = helixLayers.getByLargestContacts().getAxisAngle();
		principalRotationVector = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
	}

	/**
	 * Returns a vector perpendicular to the principal rotation vector
	 * for the alignment of structures in the xy-plane
	 * @return reference vector
	 */
	private void calcReferenceVector() {
		referenceVector = getReferenceAxisCylic();

		if (referenceVector == null) {
			logger.warn("no reference vector found. Using y-axis.");
			referenceVector = new Vector3d(Y_AXIS);
		}
		// make sure reference vector is perpendicular principal roation vector
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
		if (Math.abs(dot) < 0.00001) {
			logger.info("HelixAxisAligner: reference axis parallel");
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
