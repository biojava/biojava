package org.biojava.nbio.structure.geometry;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

/**
 * SuperPositions is a Class that provides static helper methods and an easy
 * access to the whole family of {@link SuperPosition} algorithms.
 * <p>
 * It defines a static SuperPosition object and uses it for calculation.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class SuperPositions {

	private static SuperPositionAbstract superposer = new SuperPositionQuat(
			false);
	
	/** Prevent instantiation */
	private SuperPositions(){}

	/**
	 * Use the {@link SuperPosition#superpose(Point3d[], Point3d[])} method of
	 * the default static SuperPosition algorithm contained in this Class.
	 */
	public static Matrix4d superpose(Point3d[] fixed, Point3d[] moved) {
		superposer.setCentered(false);
		return superposer.superpose(fixed, moved);
	}

	/**
	 * Use the {@link SuperPosition#superpose(Point3d[], Point3d[])} method of
	 * the default static SuperPosition algorithm contained in this Class,
	 * assuming that the point arrays are centered at the origin.
	 */
	public static Matrix4d superposeAtOrigin(Point3d[] fixed, Point3d[] moved) {
		superposer.setCentered(true);
		return superposer.superpose(fixed, moved);
	}

	/**
	 * Use the {@link SuperPosition#superposeAndTransform(Point3d[], Point3d[])}
	 * method of the default static SuperPosition algorithm contained in this
	 * Class.
	 */
	public static Matrix4d superposeAndTransform(Point3d[] fixed,
			Point3d[] moved) {
		superposer.setCentered(false);
		return superposer.superposeAndTransform(fixed, moved);
	}

	/**
	 * Use the {@link SuperPosition#superposeAndTransform(Point3d[], Point3d[])}
	 * method of the default static SuperPosition algorithm contained in this
	 * Class, assuming that the point arrays are centered at the origin.
	 */
	public static Matrix4d superposeAndTransformAtOrigin(Point3d[] fixed,
			Point3d[] moved) {
		superposer.setCentered(true);
		return superposer.superposeAndTransform(fixed, moved);
	}

	/**
	 * Use the {@link SuperPosition#getRmsd(Point3d[], Point3d[])} method of the
	 * default static SuperPosition algorithm contained in this Class.
	 */
	public static double getRmsd(Point3d[] fixed, Point3d[] moved) {
		superposer.setCentered(false);
		return superposer.getRmsd(fixed, moved);
	}

	/**
	 * Use the {@link SuperPosition#getRmsd(Point3d[], Point3d[])} method of the
	 * default static SuperPosition algorithm contained in this Class, assuming
	 * that the point arrays are centered at the origin.
	 */
	public static double getRmsdAtOrigin(Point3d[] fixed, Point3d[] moved) {
		superposer.setCentered(true);
		return superposer.getRmsd(fixed, moved);
	}
	
	public static void setDefaultSuperPosition(SuperPositionAbstract defaultAlgorithm) {
		superposer = defaultAlgorithm;
	}
}
