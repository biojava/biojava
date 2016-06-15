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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


/**
 * A cartesian product between two lists A and B is the set of all ordered pairs
 * of the elements of both sets.
 *
 * See http://en.wikipedia.org/wiki/Cartesian_product for more details.
 *
 * Since the order of the elements matters; Lists are used instead of Sets.
 *
 * Example:
 *  A = {1, 2, 3}
 *  B = {4, 5}
 *  The ordered pairs are {1, 4}, {1, 5}, {2, 4}, .., {3, 5}
 *
 * @author Peter Rose
 *
 */
public class CartesianProduct<T> {

	private List<T> list1 = Collections.emptyList();
	private List<T> list2 = Collections.emptyList();

	/**
	 * Class constructor specifying the two lists of a cartesian product.
	 */
	public CartesianProduct(List<T> list1, List<T> list2) {
		this.list1 = list1;
		this.list2 = list2;
	}

	/**
	 * Generates the list of ordered pair between two sets.
	 *
	 * @return the list of ordered pairs
	 */
	public List<OrderedPair<T>> getOrderedPairs() {
		List<OrderedPair<T>> pairs = new ArrayList<OrderedPair<T>>(list1.size()*list2.size());

		for (T element1: list1) {
			for (T element2: list2) {
				pairs.add(new OrderedPair<T>(element1, element2));
			}
		}

		return pairs;
	}
}
