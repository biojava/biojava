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
package org.biojava.nbio.structure.contact;

import java.io.Serializable;

/**
 * A Pair of objects. Implements an equals and hashCode so
 * that it can be used as key in hashes.
 * Based on the JUNG graph library implementation.
 *
 * @author duarte_j
 *
 * @param <T>
 */
public final class Pair<T> implements Serializable {

	private static final long serialVersionUID = 1L;

	private T first;
	private T second;

	/**
	 * Creates a <code>Pair</code> from the specified elements.
	 * @param value1 the first value in the new <code>Pair</code>
	 * @param value2 the second value in the new <code>Pair</code>
	 * @throws IllegalArgumentException if either argument is null
	 */
	public Pair(T value1, T value2) {
		if(value1 == null || value2 == null)
			throw new IllegalArgumentException("Pair cannot contain null values");
		first = value1;
		second = value2;
	}

	public T getFirst() {
		return first;
	}

	public T getSecond() {
		return second;
	}



	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((first == null) ? 0 : first.hashCode());
		result = prime * result + ((second == null) ? 0 : second.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		@SuppressWarnings("unchecked")
		Pair<T> other = (Pair<T>) obj;
		if (first == null) {
			if (other.first != null)
				return false;
		} else if (!first.equals(other.first))
			return false;
		if (second == null) {
			if (other.second != null)
				return false;
		} else if (!second.equals(other.second))
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "<" + first.toString() + ", " + second.toString() + ">";
	}
}
