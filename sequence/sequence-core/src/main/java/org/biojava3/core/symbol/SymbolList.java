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
package org.biojava3.core.symbol;

import java.io.Serializable;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.ListIterator;

/**
 * A symbol list is exactly the same as a normal list set up to contain symbols.
 * It is powered by a replaceable backing list which defaults to an instance of
 * {@link ArrayList}. It is possible to substitute the backing list for
 * something more efficient or better suited to the task, e.g. an interface
 * directly to disk or a database. In addition to the standard position-based
 * methods provided as part of the {@link List} interface, each method has a
 * _bio equivalent which accepts biological coordinates (1-indexed) as opposed
 * to standard coordinates (0-indexed).
 * 
 * @author Richard Holland
 * @since 3.0
 */
public class SymbolList extends AbstractList<Symbol> implements Serializable {
	
	private static final long serialVersionUID = 1L;

	private List<Symbol> backingList = new ArrayList<Symbol>();

	/**
	 * Construct a new symbol list backed by an {@link ArrayList}.
	 */
	public SymbolList() {
	}

	/**
	 * Construct a new symbol list and fill it with the symbols from the
	 * specified list. This does not wrap the specified list and operations on
	 * the constructed object do not affect the specified list.
	 * 
	 * @param symList
	 *            the symbol list to take symbols from to populate this symbol
	 *            list.
	 */
	public SymbolList(List<Symbol> symList) {
		this.addAll(symList);
	}

	/**
	 * Set the backing list for this symbol list. It will wipe out all data
	 * already in the new backing list and replace it with data already existing
	 * in this symbol list. All operations on this symbol list will directly
	 * affect the backing list, and vice versa.
	 * 
	 * @param backingList
	 *            the backing list to use.
	 */
	public void setBackingList(List<Symbol> backingList) {
		if (backingList == null) {
			throw new NullPointerException("The backing list cannot be null.");
		}
		backingList.clear();
		backingList.addAll(this);
		this.backingList = backingList;
	}

	/**
	 * Obtain the backing list currently in use. Any modifications to it will
	 * directly modify the symbol list too.
	 * 
	 * @return the backing list.
	 */
	public List<Symbol> getBackingList() {
		return this.backingList;
	}

	/**
	 * A 1-indexed equivalent of {@link add(int, Symbol)}.
	 * 
	 * @param index
	 *            the 1-indexed position.
	 * @param element
	 *            the symbol to add.
	 */
	public void add_bio(int index, Symbol element) {
		this.add(index - 1, element);
	}

	/**
	 * A 1-indexed equivalent of {@link addAll(int, Collection<? extends
	 * Symbol>)}.
	 * 
	 * @param index
	 *            the 1-indexed position.
	 * @param c
	 *            the symbols to add.
	 */
	public boolean addAll_bio(int index, Collection<? extends Symbol> c) {
		return this.addAll(index - 1, c);
	}

	/**
	 * A 1-indexed equivalent of {@link get(int)}.
	 * 
	 * @param index
	 *            the 1-indexed position.
	 * @return the symbol.
	 */
	public Symbol get_bio(int index) {
		return this.get(index - 1);
	}

	/**
	 * A 1-indexed equivalent of {@link indexOf(Object)}.
	 * 
	 * @param o
	 *            the object to look for.
	 * @return index the 1-indexed position.
	 */
	public int indexOf_bio(Object o) {
		return this.indexOf(o) + 1;
	}

	/**
	 * A 1-indexed equivalent of {@link lastIndexOf(Object)}.
	 * 
	 * @param o
	 *            the object to look for.
	 * @return index the last 1-indexed position.
	 */
	public int lastIndexOf_bio(Object o) {
		return this.lastIndexOf(o) + 1;
	}

	/**
	 * A 1-indexed equivalent of {@link remove(int)}.
	 * 
	 * @param index
	 *            the 1-indexed position.
	 * @return the symbol removed.
	 */
	public Symbol remove_bio(int index) {
		return this.remove(index - 1);
	}

	/**
	 * A 1-indexed equivalent of {@link set(int,Symbol)}.
	 * 
	 * @param index
	 *            the 1-indexed position.
	 * @param element
	 *            the symbol to add.
	 * @return the symbol replaced, if any.
	 */
	public Symbol set_bio(int index, Symbol element) {
		return this.set(index - 1, element);
	}

	public void add(int index, Symbol element) {
		this.backingList.add(index, element);
	}

	public Symbol remove(int index) {
		return this.backingList.remove(index);
	}

	public Symbol set(int index, Symbol element) {
		return this.backingList.set(index, element);
	}

	public Symbol get(int index) {
		return this.backingList.get(index);
	}

	public int size() {
		return this.backingList.size();
	}

	public ListIterator<Symbol> listIterator() {
		return this.listIterator(0);
	}

	public ListIterator<Symbol> listIterator(int index) {
		return this.listIterator(index, 0);
	}

	/**
	 * Provides a list iterator that uses 1-indexed coordinates. Otherwise
	 * identical to {@link #listIterator()}.
	 * 
	 * @return a list iterator.
	 */
	public ListIterator<Symbol> listIterator_bio() {
		return this.listIterator_bio(1);
	}

	/**
	 * Provides a list iterator that uses 1-indexed coordinates. Otherwise
	 * identical to {@link #listIterator(int)}.
	 * 
	 * @param index
	 *            the 1-indexed starting position for the iterator.
	 * @return a list iterator.
	 */
	public ListIterator<Symbol> listIterator_bio(int index) {
		return this.listIterator(index, 1);
	}

	/**
	 * Provides a list iterator that uses offest-indexed coordinates.
	 * 
	 * @param index
	 *            the offset-indexed starting position for the iterator.
	 * @param offset
	 *            the offset for indexing (usually 0 or 1).
	 * @return a list iterator.
	 */
	private ListIterator<Symbol> listIterator(final int index, final int offset) {
		final ListIterator<Symbol> def = super.listIterator(index - offset);
		return new ListIterator<Symbol>() {

			public void add(Symbol e) {
				def.add(e);
			}

			public boolean hasNext() {
				return def.hasNext();
			}

			public boolean hasPrevious() {
				return def.hasPrevious();
			}

			public Symbol next() {
				return def.next();
			}

			public int nextIndex() {
				return def.nextIndex() + offset;
			}

			public Symbol previous() {
				return def.previous();
			}

			public int previousIndex() {
				return def.previousIndex() + offset;
			}

			public void remove() {
				def.remove();
			}

			public void set(Symbol e) {
				def.set(e);
			}
		};
	}
}
