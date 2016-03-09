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
package org.biojava.nbio.structure.quaternary;

/**
 * An ordered pair represents a component of a cartesian product. The cartesian product
 * of two sets A and B is the set of all ordered pairs of the elements of both sets.
 *
 * See http://en.wikipedia.org/wiki/Cartesian_product for more details.
 *
 * Example:
 *  A = {1, 2, 3}
 *  B = {4, 5}
 *  The ordered pairs are {1, 4}, {1, 5}, {2, 4}, .., {3, 5}
 *
 * @author Peter Rose
 *
 * @param <T>
 */
public class OrderedPair<T> {
	T element1;
	T element2;

	/**
	 * Class constructor specifying the two elements of an ordered pair.
	 */
	OrderedPair(T element1, T element2) {
		this.element1 = element1;
		this.element2 = element2;
	}

	/**
	 * @return element1 the first element of an ordered pair
	 */
	public T getElement1() {
		return element1;
	}

	/**
	 * Sets the first element of an ordered pair.
	 *
	 * @param element1 the first element of an ordered pair
	 */
	public void setElement1(T element1) {
		this.element1 = element1;
	}

	/**
	 * @return element2 the second element of an ordered pair
	 */
	public T getElement2() {
		return element2;
	}

	/**
	 * Sets the second element of an ordered pair.
	 *
	 * @param element2 the second element of an ordered pair
	 */
	public void setElement2(T element2) {
		this.element2 = element2;
	}

	@Override
	public String toString() {
		return "[" + element1.toString() + "," +
		element2.toString() + "]";
	}
}
