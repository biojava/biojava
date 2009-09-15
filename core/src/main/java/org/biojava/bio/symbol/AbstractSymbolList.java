/*
 * BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 * http://www.biojava.org
 *
 */

package org.biojava.bio.symbol;

import java.io.Serializable;
import java.util.AbstractList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <p>
 * Abstract helper implementation of the SymbolList core interface. To produce a
 * concrete SymbolList implementation, you need only implement the
 * <code>getAlphabet</code>, <code>length</code> and <code>symbolAt</code>
 * methods. Iterators and sublists are handled for you automatically.
 * </p>
 * 
 * <p>
 * This class makes many custom SymbolList implementations very quick to
 * implement. See org.biojava.bio.seq.tools.ComplementSymbolList for an example
 * of this.
 * </p>
 * 
 * <p>
 * To make a mutable SymbolList, override the implementation of edit to perform
 * the apropreate edit. If your implementation of SymbolList is a view onto an
 * underlying SymbolList, then you must forward any apropreate edit requests to
 * that list, and forward all events from the underlying list to your listeners.
 * </p>
 * 
 * @author Thomas Down
 * @author Matthew Pocock
 */

public abstract class AbstractSymbolList extends AbstractChangeable implements
		SymbolList {
	protected AbstractSymbolList() {
	}

	public Iterator<Symbol> iterator() {
		return new SymbolIterator(1, length());
	}

	public SymbolList subList(int start, int end) {
		if (start < 1 || end > length()) {
			throw new IndexOutOfBoundsException("Sublist index out of bounds "
					+ length() + ":" + start + "," + end);
		}

		if (end < start) {
			throw new IllegalArgumentException(
					"end must not be lower than start: start=" + start
							+ ", end=" + end);
		}
		return new SubList(start, end);
	}

	public List toList() {
		return new ListView(this);
	}

	public String seqString() {
		try {
			SymbolTokenization toke = getAlphabet().getTokenization("default");
			return toke.tokenizeSymbolList(this);
		} catch (BioException ex) {
			throw new BioRuntimeException("Couldn't tokenize sequence", ex);
		}
	}

	public String subStr(int start, int end) {
		return subList(start, end).seqString();
	}

	public void edit(Edit edit) throws IllegalAlphabetException,
			ChangeVetoException {
		throw new ChangeVetoException("AbstractSymbolList is immutable");
	}

	/**
	 * Provides logical equality for two SymbolLists that share the same list of
	 * canonical symbols
	 */
	public boolean equals(Object o) {
		if (this == o)
			return true;// just for optimality
		if (o == null)
			return false;
		if (!(o instanceof SymbolList))
			return false;
		return compare(this, (SymbolList) o);
	}

	private volatile int hashCode = 0; // for caching purposes

	public int hashCode() {
		if (hashCode == 0) {
			int result = 17;
			for (Iterator i = iterator(); i.hasNext();) {
				result = 37 * result + i.next().hashCode();
			}
			hashCode = result;
		}
		return hashCode;
	}

	public String toString() {
		return super.toString() + " length: " + length();
	}

	private boolean compare(SymbolList sl1, SymbolList sl2) {
		if (!sl1.getAlphabet().equals(sl2.getAlphabet())) {
			return false;
		}

		if (sl1.length() != sl2.length()) {
			return false;
		}

		Iterator si1 = sl1.iterator();
		Iterator si2 = sl2.iterator();
		while (si1.hasNext()) {
			if (!(si1.next() == si2.next())) {
				return false;
			}
		}
		return true;
	}

	/**
	 * <p>
	 * An Iterator over each Symbol in a range of a SymbolList.
	 * </p>
	 * 
	 * <p>
	 * Objects of this type are returned by
	 * <code>AbstractSymbolList.iterator</code>.
	 * </p>
	 * 
	 * @author Thomas Down
	 */

	private class SymbolIterator implements Iterator<Symbol>, Serializable {
		/**
		 * Generated Serial Version ID.
		 */
		private static final long serialVersionUID = -8542733313998842964L;
		private int max;
		private int pos;

		/**
		 * Construct a SymbolIterator object that will return the symbols from
		 * min to max inclusive.
		 * 
		 * @param min
		 *            the first index to return
		 * @param max
		 *            the last index to return
		 */

		public SymbolIterator(int min, int max) {
			this.max = max;
			pos = min;
		}

		protected SymbolIterator() {
		}

		public boolean hasNext() {
			return (pos <= max);
		}

		public Symbol next() {
			if (pos > max) {
				throw new NoSuchElementException();
			}
			return symbolAt(pos++);
		}

		public void remove() {
			throw new UnsupportedOperationException();
		}
	}

	/**
	 * <p>
	 * Implements a list view of a SymbolList.
	 * </p>
	 * 
	 * <p>
	 * Objects of this type are instantiated by
	 * <code>AbstractSymbolList.subList</code>.
	 * </p>
	 * 
	 * @author Thomas Down
	 * @author Matthew Pocock
	 */

	private class SubList extends AbstractChangeable implements SymbolList,
			Serializable {
		private int start, end;

		private transient EditTranslater editTranslater = null;
		private transient ChangeForwarder annotationForwarder = null;

		public SubList(int start, int end) {
			this.start = start;
			this.end = end;
		}

		public Alphabet getAlphabet() {
			return AbstractSymbolList.this.getAlphabet();
		}

		public Iterator iterator() {
			return new SymbolIterator(start, end);
		}

		public int length() {
			return end - start + 1;
		}

		public Symbol symbolAt(int pos) {
			if (pos < 1 || pos > length()) {
				throw new IndexOutOfBoundsException(
						"Symbol index out of bounds " + length() + ":" + pos);
			}
			return AbstractSymbolList.this.symbolAt(pos + start - 1);
		}

		public SymbolList subList(int sstart, int send) {
			if (sstart < 1 || send > length()) {
				throw new IndexOutOfBoundsException(
						"Sublist index out of bounds " + length() + ":"
								+ sstart + "," + send);
			}

			if (send < sstart) {
				throw new IndexOutOfBoundsException(
						"Requested end must not be lower than start: start="
								+ sstart + ", end=" + send);
			}

			return new SubList(sstart + start - 1, send + start - 1);
		}

		public String seqString() {
			try {
				SymbolTokenization toke = getAlphabet()
						.getTokenization("token");
				return toke.tokenizeSymbolList(this);
			} catch (BioException ex) {
				throw new BioRuntimeException("Couldn't tokenize sequence", ex);
			}
		}

		public String subStr(int start, int end) {
			return subList(start, end).seqString();
		}

		public List toList() {
			return new ListView(this);
		}

		// fixme: doesn't do range checking on edit object
		public void edit(Edit edit) throws IllegalAlphabetException,
				ChangeVetoException {
			AbstractSymbolList.this.edit(new Edit(edit.pos + this.start - 1,
					edit.length, edit.replacement));
		}

		protected ChangeSupport getChangeSupport(ChangeType changeType) {
			ChangeSupport cs = super.getChangeSupport(changeType);

			if ((SymbolList.EDIT.isMatchingType(changeType) || changeType
					.isMatchingType(SymbolList.EDIT))
					&& (editTranslater == null)) {
				editTranslater = new EditTranslater(this, cs, start, end);
				AbstractSymbolList.this.addChangeListener(editTranslater,
						SymbolList.EDIT);
			}

			if (((changeType == null) || (changeType == Annotation.PROPERTY))
					&& (annotationForwarder == null)) {
				annotationForwarder = new ChangeForwarder.Retyper(this, cs,
						Annotation.PROPERTY);
				AbstractSymbolList.this.addChangeListener(annotationForwarder,
						Annotation.PROPERTY);
			}

			return cs;
		}

		/**
		 * Provides logical equality for two SymbolLists that share the same
		 * list of canonical symbols
		 */
		public boolean equals(Object o) {
			if (this == o)
				return true;// just for optimality
			if (o == null)
				return false;
			if (!(o instanceof SymbolList))
				return false;
			return compare(this, (SymbolList) o);
		}

		private volatile int hashCode = 0; // for caching purposes

		public int hashCode() {
			if (hashCode == 0) {
				int result = 17;
				for (Iterator i = iterator(); i.hasNext();) {
					result = 37 * result + i.next().hashCode();
				}
				hashCode = result;
			}
			return hashCode;
		}

		public String toString() {
			return super.toString() + " start: " + start + " end: " + end /*
																		 * (+
																		 * " parent: "
																		 * +
																		 * AbstractSymbolList
																		 * .
																		 * this.
																		 * toString
																		 * ()
																		 * [causes
																		 * infinite
																		 * loop
																		 * --
																		 * odd]
																		 */;
		}
	}

	/**
	 * <p>
	 * Implements a list view of a SymbolList.
	 * </p>
	 * 
	 * <p>
	 * Objects of this type are instantiated by
	 * <code>AbstractSymbolList.asList</code>.
	 * </p>
	 * 
	 * @author Thomas Down
	 */

	private static class ListView extends AbstractList implements Serializable {
		private SymbolList rl;

		/**
		 * Build a new ListView to view a SymbolList as a list.
		 * 
		 * @param symList
		 *            the SymbolList to view
		 */

		ListView(SymbolList symList) {
			this.rl = symList;
		}

		public Object get(int pos) {
			return rl.symbolAt(pos + 1);
		}

		public int size() {
			return rl.length();
		}
	}

	/**
	 * This adapter screens all edit events to see if they overlap with a window
	 * of interest. If they do, then a new edit event is built for the
	 * overlapping region and pased on to all listeners.
	 * 
	 * @author Matthew Pocock
	 */
	public class EditScreener extends ChangeForwarder {
		protected final int min;
		protected final int max;

		public EditScreener(Object source, ChangeSupport cs, int min, int max) {
			super(source, cs);
			this.min = min;
			this.max = max;
		}

		protected ChangeEvent generateEvent(ChangeEvent ce) {
			ChangeType ct = ce.getType();
			if (ct == SymbolList.EDIT) {
				Object change = ce.getChange();
				if ((change != null) && (change instanceof Edit)) {
					Edit edit = (Edit) change;
					int start = edit.pos;
					int end = start + edit.length - 1; // inclusive
					if ((start <= max) && (end >= min)) {
						// they overlap
						return new ChangeEvent(getSource(), ct, edit, null, ce);
					}
				}
			}
			return null;
		}
	}

	/**
	 * This translates edit events that fall within a window into window
	 * co-ordinates.
	 * 
	 * @author Matthew Pocock
	 */
	public class EditTranslater extends EditScreener {
		public EditTranslater(Object source, ChangeSupport cs, int min, int max) {
			super(source, cs, min, max);
		}

		protected ChangeEvent generateEvent(ChangeEvent ce) {
			ce = super.generateEvent(ce);
			if (ce != null) {
				Edit edit = (Edit) ce.getChange();
				return new ChangeEvent(ce.getSource(), ce.getType(), new Edit(
						edit.pos - min, edit.length, edit.replacement), null,
						ce.getChainedEvent());
			}
			return null;
		}
	}
}
