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
package org.biojava.nbio.ontology.utils;

import java.io.Serializable;
import java.util.*;

/**
 * Lightweight implementation of Map which uses little memory to store a
 * small number of mappings, at the expense of scalability.  Not recommended
 * for more than 20-30 mappings.
 *
 * <p>
 * This implementation has the useful property that the iteration order is
 * the same as the order in which mappings are added.
 * </p>
 *
 * @author Thomas Down
 * @author Len Trigg
 */

public class SmallMap extends AbstractMap implements Serializable {
	private Object[] mappings = null;
	private int numMappings = 0;

	public SmallMap() {
		super();
	}

	public SmallMap(int size) {
		super();
		mappings = new Object[size * 2];
	}

	public SmallMap(Map m) {
		this(m.size());
		for (Iterator i = m.entrySet().iterator(); i.hasNext(); ) {
			Map.Entry me = (Map.Entry) i.next();
			put(me.getKey(), me.getValue());
		}
	}

	/**
	 * @throws NullPointerException if key is null
	 */
	@Override
	public Object get(Object key) {
		// Doesn't actually need to check if mappings is null, since numMappings
		// will necessarily be zero.

		int keyHash = key.hashCode();
		for (int i = 0; i < numMappings * 2; i += 2) {
			if (keyHash == mappings[i].hashCode() && key.equals(mappings[i])) {
				return mappings[i + 1];
			}
		}
		return null;
	}

	/**
	 * @throws NullPointerException if key is null
	 */
	@Override
	public Object put(Object key, Object value) {
		int keyHash = key.hashCode();

		for (int i = 0; i < numMappings * 2; i += 2) {
			if (keyHash == mappings[i].hashCode() && key.equals(mappings[i])) {
				Object oldValue = mappings[i+1];
				mappings[i+1] = value;
				return oldValue;
			}
		}

		int newPos = numMappings * 2;
		int oldLength = 0;
		if (mappings != null) {
			oldLength = mappings.length;
		}
		if (newPos + 1 >= oldLength) {
			Object[] newMappings = new Object[oldLength + 6];
			if (oldLength > 0) {
				System.arraycopy(mappings, 0, newMappings, 0, mappings.length);
			}
			mappings = newMappings;
		}

		mappings[newPos] = key;
		mappings[newPos + 1] = value;
		numMappings++;

		return null;
	}

	@Override
	public Set keySet() {
		return new KeySet();
	}

	@Override
	public Set entrySet() {
		return new EntrySet();
	}

	@Override
	public int size() {
		return numMappings;
	}

	@Override
	public boolean containsKey(Object key) {
		int keyHash = key.hashCode();
		for (int i = 0; i < numMappings * 2; i += 2) {
			if (keyHash == mappings[i].hashCode() && key.equals(mappings[i])) {
				return true;
			}
		}
		return false;
	}


	// num ranges from 1 to numMappings
	private void removeMapping(int num) {
		if (num < numMappings) {
			System.arraycopy(mappings, num * 2, mappings, (num - 1) * 2, (numMappings - num) * 2);
		}
		mappings[numMappings * 2 - 1] = null;
		mappings[numMappings * 2 - 2] = null;
		numMappings--;
	}

	private class KeySet extends AbstractSet {
		@Override
		public Iterator iterator() {
			return new KeyIterator();
		}

		@Override
		public int size() {
			return numMappings;
		}
	}

	private class KeyIterator implements Iterator {
		int pos = 0;

		@Override
		public boolean hasNext() {
			return pos < numMappings;
		}

		@Override
		public Object next() {
			if (pos >= numMappings) {
				throw new NoSuchElementException();
			}
			int offset = pos * 2;
			++pos;
			return mappings[offset];
		}

		@Override
		public void remove() {
			removeMapping(pos);
			pos -= 1;
		}
	}

	private class EntrySet extends AbstractSet {
		@Override
		public Iterator iterator() {
			return new EntryIterator();
		}

		@Override
		public int size() {
			return numMappings;
		}
	}

	private class EntryIterator implements Iterator {
		int pos = 0;

		@Override
		public boolean hasNext() {
			return pos < numMappings;
		}

		@Override
		public Object next() {
			if (pos >= numMappings) {
				throw new NoSuchElementException();
			}
			int offset = pos * 2;
			++pos;
			return new MapEntry(offset);
		}

		@Override
		public void remove() {
			removeMapping(pos);
			pos -= 1;
		}
	}

	private class MapEntry implements Map.Entry {
		private int offset;

		private MapEntry(int offset) {
			this.offset = offset;
		}

		@Override
		public Object getKey() {
			return mappings[offset];
		}

		@Override
		public Object getValue() {
			return mappings[offset + 1];
		}

		@Override
		public Object setValue(Object v) {
			Object oldValue = mappings[offset + 1];
			mappings[offset + 1] = v;
			return oldValue;
		}

		@Override
		public boolean equals(Object o) {
			if (! (o instanceof Map.Entry)) {
				return false;
			}

			Map.Entry mo = (Map.Entry) o;
			return ((getKey() == null ? mo.getKey() == null : getKey().equals(mo.getKey())) &&
			(getValue() == null ? mo.getValue() == null : getValue().equals(mo.getValue())));
		}

		@Override
		public int hashCode() {
			return (getKey() == null ? 0 : getKey().hashCode()) ^ (getValue() == null ? 0 : getValue().hashCode());
		}
	}
}
