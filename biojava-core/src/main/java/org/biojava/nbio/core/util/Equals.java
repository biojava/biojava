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
package org.biojava.nbio.core.util;

/**
 * A set of helper methods which return true if the two parameters are
 * equal to each other.
 *
 * @see #classEqual(Object, Object) For how to use this class
 *
 * @author ayates
 */
public class Equals {

	public static boolean equal(int one, int two) {
		return one == two;
	}

	public static boolean equal(long one, long two) {
		return (one == two);
	}

	public static boolean equal(boolean one, boolean two) {
		return one == two;
	}

	/**
	 * Does not compare class types.
	 * @see #classEqual(Object, Object)
	 */
	public static boolean equal(Object one, Object two) {
		return one == null && two == null || !(one == null || two == null) && (one == two || one.equals(two));
	}

	/**
	 * This method should be called before beginning any equals methods. In order
	 * to return true the method:
	 *
	 * <ol>
	 * <li>The two given objects are the same instance using ==. This also means
	 * if both Objects are null then this method will return true (well
	 * technically they are equal)</li>
	 * <li>Tests that neither object is null</li>
	 * <li>The the two classes from the objects are equal using ==</li>
	 * </ol>
	 *
	 * The boilerplate using this method then becomes:
	 *
	 * <pre>
	 * boolean equals = false;
	 * if (EqualsHelper.classEqual(this, obj)) {
	 *   TargetClass casted = (TargetClass) obj;
	 *   equals = (EqualsHelper.equal(this.getId(), casted.getId()) &amp;&amp; EqualsHelper
	 *       .equal(this.getName(), casted.getName()));
	 * }
	 * return equals;
	 * </pre>
	 *
	 * @param one
	 *          The first object to test
	 * @param two
	 *          The second object to test
	 * @return A boolean indicating if the logic agrees that these two objects are
	 *         equal at the class level
	 */
	public static boolean classEqual(Object one, Object two) {
		return one == two || !(one == null || two == null) && one.getClass() == two.getClass();
	}
}
