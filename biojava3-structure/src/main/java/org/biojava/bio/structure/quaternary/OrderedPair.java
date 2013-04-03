package org.biojava.bio.structure.quaternary;

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
	
	public String toString() {
		return "[" + element1.toString() + "," + 
		element2.toString() + "]";
	}
}
