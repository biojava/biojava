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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.NoSuchElementException;

/**
 * A tagged symbol list behaves identically to a normal {@link SymbolList},
 * except that it maintains in parallel a set of tags of equal size. By default
 * every tag is null. It contains methods for changing and interrogating the tag
 * at each position of the symbol list, and for iterating over the set of tags.
 * It is also possible to iterate over a set of {@link TaggedSymbol}
 * representations containing both the symbol and the tag from each position.
 * Writing a tagged symbol list out using the standard {@link SymbolListFormat}
 * methods does not preserve the tags.
 * 
 * @author Richard Holland
 * @since 3.0
 * @param T
 *            the type of object to use for the tags.
 */
public class TaggedSymbolList<T> extends SymbolList implements Serializable {

	private static final long serialVersionUID = 1L;
	
	private List<T> tags = new ArrayList<T>();

	/**
	 * Construct a new, empty, tagged symbol list.
	 */
	public TaggedSymbolList() {
	}

	/**
	 * Construct a new tagged symbol list and fill it with the symbols from the
	 * specified list, with null tags at every position. This does not wrap the
	 * specified list and operations on the constructed object do not affect the
	 * specified list.
	 * 
	 * @param symList
	 *            the symbol list to take symbols from to populate the tagged
	 *            symbol list.
	 */
	public TaggedSymbolList(List<Symbol> symList) {
		this.addAll(symList);
	}

	@Override
	public void add(int index, Symbol element) {
		this.tags.add(index, null);
		super.add(index, element);
	}

	@Override
	public Symbol remove(int index) {
		this.tags.remove(index);
		return super.remove(index);
	}

	/**
	 * Adds a symbol and associates a tag with it.
	 * 
	 * @param e
	 *            the symbol to add.
	 * @param tag
	 *            the tag to associated with it.
	 */
	public void add(Symbol e, T tag) {
		this.add(0, e, tag);
	}

	/**
	 * Adds a symbol and associates a tag with it.
	 * 
	 * @param index
	 *            the 0-indexed position to insert at.
	 * @param e
	 *            the symbol to add.
	 * @param tag
	 *            the tag to associated with it.
	 * @return {@code true} (as specified by {@link Collection#add})
	 */
	public void add(int index, Symbol e, T tag) {
		this.tags.add(index, tag);
		this.add(index, e);
	}

	/**
	 * Adds a symbol and associates a tag with it.
	 * 
	 * @param index
	 *            the 1-indexed position to insert at.
	 * @param e
	 *            the symbol to add.
	 * @param tag
	 *            the tag to associated with it.
	 * @return {@code true} (as specified by {@link Collection#add})
	 */
	public void add_bio(int index, Symbol e, T tag) {
		this.add(index - 1, e, tag);
	}

	/**
	 * Adds symbols and associates the same tag with each one.
	 * 
	 * @param c
	 *            the symbols to add.
	 * @param tag
	 *            the tag to associate with the added symbols.
	 * @return {@code true} (as specified by {@link Collection#addAll})
	 */
	public boolean addAll(Collection<? extends Symbol> c, T tag) {
		return this.addAll(0, c, tag);
	}

	/**
	 * Adds symbols and associates each one with the corresponding tag. The two
	 * collections must be the same size and correspondence is determined by
	 * iteration.
	 * 
	 * @param c
	 *            the symbols to add.
	 * @param tags
	 *            the tags to associate with the added symbols.
	 * @return {@code true} (as specified by {@link Collection#addAll})
	 */
	public boolean addAll(Collection<? extends Symbol> c, Collection<T> tags) {
		return this.addAll(0, c, tags);
	}

	/**
	 * Adds symbols and associates the same tag with each one.
	 * 
	 * @param index
	 *            the 0-indexed position to insert at.
	 * @param c
	 *            the symbols to add.
	 * @param tag
	 *            the tag to associate with the added symbols.
	 * @return {@code true} (as specified by {@link Collection#addAll})
	 */
	public boolean addAll(int index, Collection<? extends Symbol> c, T tag) {
		for (int i = 0; i < c.size(); i++) {
			this.tags.add(index + i, tag);
		}
		return super.addAll(index, c);
	}

	/**
	 * Adds symbols and associates each one with the corresponding tag. The two
	 * collections must be the same size and correspondence is determined by
	 * iteration.
	 * 
	 * @param index
	 *            the 0-indexed position to insert at.
	 * @param c
	 *            the symbols to add.
	 * @param tags
	 *            the tags to associate with the added symbols.
	 * @return {@code true} (as specified by {@link Collection#addAll})
	 */
	public boolean addAll(int index, Collection<? extends Symbol> c,
			Collection<T> tags) {
		if (tags.size() != c.size()) {
			throw new IndexOutOfBoundsException(
					"Tag and symbol collections must be the same size.");
		}
		this.tags.addAll(index, tags);
		return super.addAll(index, c);
	}

	/**
	 * Adds symbols and associates the same tag with each one.
	 * 
	 * @param index
	 *            the 1-indexed position to insert at.
	 * @param c
	 *            the symbols to add.
	 * @param tag
	 *            the tag to associate with the added symbols.
	 * @return {@code true} (as specified by {@link Collection#addAll})
	 */
	public boolean addAll_bio(int index, Collection<? extends Symbol> c, T tag) {
		return this.addAll(index - 1, c, tag);
	}

	/**
	 * Adds symbols and associates each one with the corresponding tag. The two
	 * collections must be the same size and correspondence is determined by
	 * iteration.
	 * 
	 * @param index
	 *            the 1-indexed position to insert at.
	 * @param c
	 *            the symbols to add.
	 * @param tags
	 *            the tags to associate with the added symbols.
	 * @return {@code true} (as specified by {@link Collection#addAll})
	 */
	public boolean addAll_bio(int index, Collection<? extends Symbol> c,
			Collection<T> tags) {
		return this.addAll(index - 1, c, tags);
	}

	/**
	 * Retrieves the tag at the specified position.
	 * 
	 * @param index
	 *            the 0-indexed position.
	 * @return the tag.
	 */
	public T getTag(int index) {
		return this.tags.get(index);
	}

	/**
	 * Retrieves the tag at the specified position.
	 * 
	 * @param index
	 *            the 1-indexed position.
	 * @return the tag.
	 */
	public T getTag_bio(int index) {
		return this.getTag(index - 1);
	}

	/**
	 * Sets the symbol and tag at the specified position.
	 * 
	 * @param index
	 *            the 0-indexed position.
	 * @param element
	 *            the symbol.
	 * @param tag
	 *            the tag.
	 * @return the element previously at the specified position.
	 */
	public Symbol set(int index, Symbol element, T tag) {
		this.tags.set(index, tag);
		return super.set(index, element);
	}

	/**
	 * Sets the symbol and tag at the specified position.
	 * 
	 * @param index
	 *            the 1-indexed position.
	 * @param element
	 *            the symbol.
	 * @param tag
	 *            the tag.
	 * @return the element previously at the specified position.
	 */
	public Symbol set_bio(int index, Symbol element, T tag) {
		return this.set(index - 1, element, tag);
	}

	/**
	 * Sets the tag at the specified position.
	 * 
	 * @param index
	 *            the 0-indexed position.
	 * @param tag
	 *            the tag.
	 * @return the tag previously at the specified position.
	 */
	public T setTag(int index, T tag) {
		return this.tags.set(index, tag);
	}

	/**
	 * Sets the tag at the specified position.
	 * 
	 * @param index
	 *            the 1-indexed position.
	 * @param tag
	 *            the tag.
	 * @return the tag previously at the specified position.
	 */
	public T setTag_bio(int index, T tag) {
		return this.setTag(index - 1, tag);
	}

	/**
	 * Fills the given range with the given tag.
	 * 
	 * @param start
	 *            the 0-indexed start.
	 * @param length
	 *            how many symbols to assign the tag to.
	 * @param tag
	 *            the tag to assign.
	 */
	public void setTagRange(int start, int length, T tag) {
		for (int i = start; i < start + length; i++) {
			this.setTag(i, tag);
		}
	}

	/**
	 * Fills the given range with the given tag.
	 * 
	 * @param start
	 *            the 1-indexed start.
	 * @param length
	 *            how many symbols to assign the tag to.
	 * @param tag
	 *            the tag to assign.
	 */
	public void setTagRange_bio(int start, int length, T tag) {
		this.setTagRange(start - 1, length, tag);
	}

	/**
	 * Fills the given range with the given collection of tags, which must be
	 * equal in size to the length of the range. Tag order is determined by
	 * iteration.
	 * 
	 * @param start
	 *            the 0-indexed start.
	 * @param length
	 *            how many symbols to assign the tag to.
	 * @param tags
	 *            the tags to assign.
	 */
	public void setTagRange(int start, int length, Collection<T> tags) {
		if (tags.size() != length) {
			throw new IndexOutOfBoundsException(
					"Cannot set range of tags of different length to collection size.");
		}
		for (Iterator<T> i = tags.iterator(); i.hasNext();) {
			this.setTag(start++, i.next());
		}
	}

	/**
	 * Fills the given range with the given collection of tags, which must be
	 * equal in size to the length of the range. Tag order is determined by
	 * iteration.
	 * 
	 * @param start
	 *            the 1-indexed start.
	 * @param length
	 *            how many symbols to assign the tag to.
	 * @param tags
	 *            the tags to assign.
	 */
	public void setTagRange_bio(int start, int length, Collection<T> tags) {
		this.setTagRange(start - 1, length, tags);
	}

	/**
	 * Iterate over all the tags in order. Removing from the iterator will
	 * remove the entire tag+symbol pair from the tagged symbol list.
	 * 
	 * @return an iterator over each tag.
	 */
	public Iterator<T> tagIterator() {
		return new Iterator<T>() {

			private int pos = 0;
			private boolean removed = true;

			public boolean hasNext() {
				return this.pos < size();
			}

			public T next() {
				if (!this.hasNext()) {
					throw new NoSuchElementException();
				}
				if (!this.removed) {
					this.pos++;
				}
				this.removed = false;
				return TaggedSymbolList.this.getTag(pos);
			}

			public void remove() {
				if (this.removed) {
					throw new IllegalStateException();
				}
				TaggedSymbolList.this.remove(this.pos);
				this.removed = true;
			}
		};
	}

	/**
	 * Iterate over all the tags in order. Removing from the iterator will
	 * remove the entire tag+symbol pair from the tagged symbol list.
	 * 
	 * @return an iterator over the tags.
	 */
	public ListIterator<T> tagListIterator() {
		return this.tagListIterator(0);
	}

	/**
	 * Iterate over all the tags in order. Removing from the iterator will
	 * remove the entire tag+symbol pair from the tagged symbol list.
	 * 
	 * @param index
	 *            the 0-indexed position to start iterating from.
	 * @return an iterator over the tags.
	 */
	public ListIterator<T> tagListIterator(final int index) {
		return this.tagListIterator(index, 0);
	}

	/**
	 * Iterate over all the tags in order. Removing from the iterator will
	 * remove the entire tag+symbol pair from the tagged symbol list. Uses
	 * 1-indexed coordinates for the iterator.
	 * 
	 * @return an iterator over the tags.
	 */
	public ListIterator<T> tagListIterator_bio() {
		return this.tagListIterator_bio(1);
	}

	/**
	 * Iterate over all the tags in order. Removing from the iterator will
	 * remove the entire tag+symbol pair from the tagged symbol list.
	 * 
	 * @param index
	 *            the 1-indexed position to start iterating from.
	 * @return an iterator over the tags.
	 */
	public ListIterator<T> tagListIterator_bio(final int index) {
		return this.tagListIterator(index, 1);
	}

	/**
	 * Iterate over all the tags in order. Removing from the iterator will
	 * remove the entire tag+symbol pair from the tagged symbol list.
	 * 
	 * @param index
	 *            the offset-indexed position to start iterating from.
	 * @param offset
	 *            the offset to use for indexing, usually 0 or 1.
	 * @return an iterator over the tags.
	 */
	private ListIterator<T> tagListIterator(final int index, final int offset) {
		return new ListIterator<T>() {

			private int pos = index - offset;
			private boolean removed = true;

			public void add(T e) {
				throw new UnsupportedOperationException(
						"Cannot add to tag list iterator.");
			}

			public boolean hasPrevious() {
				return this.pos > 0;
			}

			public int nextIndex() {
				return this.pos + 1 + offset;
			}

			public T previous() {
				this.pos--;
				this.removed = false;
				return TaggedSymbolList.this.getTag(pos);
			}

			public int previousIndex() {
				return this.pos - 1 + offset;
			}

			public void set(T e) {
				TaggedSymbolList.this.setTag(this.pos, e);
			}

			public boolean hasNext() {
				return this.pos < size();
			}

			public T next() {
				if (!this.removed) {
					this.pos++;
				}
				this.removed = false;
				return TaggedSymbolList.this.getTag(pos);
			}

			public void remove() {
				if (!this.removed) {
					TaggedSymbolList.this.remove(this.pos);
					this.removed = true;
				}
			}
		};
	}

	/**
	 * Iterate over all the symbols and tags and return each pair as a
	 * {@link TaggedSymbol} instance. Modifying the iterator or the returned
	 * pairs in any way will directly modify the underlying tagged symbol list.
	 * 
	 * @return the iterator over the symbol and tag pairs.
	 */
	public Iterator<TaggedSymbol<T>> taggedSymbolIterator() {
		return new Iterator<TaggedSymbol<T>>() {

			private int pos = 0;
			private boolean removed = true;

			public boolean hasNext() {
				return this.pos < size();
			}

			public TaggedSymbol<T> next() {
				if (!this.removed) {
					this.pos++;
				}
				this.removed = false;
				return new TaggedSymbol<T>() {	
					
					private static final long serialVersionUID = 1L;
				
					private final int symPos = pos;

					public Symbol getSymbol() {
						return TaggedSymbolList.this.get(symPos);
					}

					public T getTag() {
						return TaggedSymbolList.this.getTag(symPos);
					}

					public Symbol setSymbol(Symbol sym) {
						return TaggedSymbolList.this.set(symPos, sym);
					}

					public T setTag(T tag) {
						return TaggedSymbolList.this.setTag(symPos, tag);
					}
				};
			}

			public void remove() {
				if (!this.removed) {
					TaggedSymbolList.this.remove(this.pos);
					this.removed = true;
				}
			}
		};
	}

	/**
	 * Iterate over all the symbols and tags and return each pair as a
	 * {@link TaggedSymbol} instance. Modifying the iterator or the returned
	 * pairs in any way will directly modify the underlying tagged symbol list.
	 * 
	 * @return the iterator over the symbol and tag pairs.
	 */
	public ListIterator<TaggedSymbol<T>> taggedSymbolListIterator() {
		return this.taggedSymbolListIterator(0);
	}

	/**
	 * Iterate over all the symbols and tags and return each pair as a
	 * {@link TaggedSymbol} instance. Modifying the iterator or the returned
	 * pairs in any way will directly modify the underlying tagged symbol list.
	 * 
	 * @param index
	 *            the 0-indexed position to start iterating from.
	 * @return the iterator over the symbol and tag pairs.
	 */
	public ListIterator<TaggedSymbol<T>> taggedSymbolListIterator(
			final int index) {
		return this.taggedSymbolListIterator(index, 0);
	}

	/**
	 * Iterate over all the symbols and tags and return each pair as a
	 * {@link TaggedSymbol} instance. Modifying the iterator or the returned
	 * pairs in any way will directly modify the underlying tagged symbol list.
	 * Uses 1-indexed coordinates for the iterator.
	 * 
	 * @return the iterator over the symbol and tag pairs.
	 */
	public ListIterator<TaggedSymbol<T>> taggedSymbolListIterator_bio() {
		return this.taggedSymbolListIterator_bio(1);
	}

	/**
	 * Iterate over all the symbols and tags and return each pair as a
	 * {@link TaggedSymbol} instance. Modifying the iterator or the returned
	 * pairs in any way will directly modify the underlying tagged symbol list.
	 * 
	 * @param index
	 *            the 1-indexed position to start iterating from.
	 * @return the iterator over the symbol and tag pairs.
	 */
	public ListIterator<TaggedSymbol<T>> taggedSymbolListIterator_bio(
			final int index) {
		return this.taggedSymbolListIterator(index, 1);
	}

	/**
	 * Iterate over all the symbols and tags and return each pair as a
	 * {@link TaggedSymbol} instance. Modifying the iterator or the returned
	 * pairs in any way will directly modify the underlying tagged symbol list.
	 * 
	 * @param index
	 *            the offset-indexed position to start iterating from.
	 * @param offset
	 *            the offset to use for indexing, usually 0 or 1.
	 * @return the iterator over the symbol and tag pairs.
	 */
	private ListIterator<TaggedSymbol<T>> taggedSymbolListIterator(
			final int index, final int offset) {
		return new ListIterator<TaggedSymbol<T>>() {

			private int pos = index - offset;
			private boolean removed = true;

			public void add(TaggedSymbol<T> e) {
				throw new UnsupportedOperationException(
						"Cannot add to tagged symbol list iterator.");
			}

			public boolean hasPrevious() {
				return this.pos > 0;
			}

			public int nextIndex() {
				return this.pos + 1 + offset;
			}

			public TaggedSymbol<T> previous() {
				this.pos--;
				this.removed = false;
				return new TaggedSymbol<T>() {
					
					private static final long serialVersionUID = 1L;
					
					private final int symPos = pos;

					public Symbol getSymbol() {
						return TaggedSymbolList.this.get(symPos);
					}

					public T getTag() {
						return TaggedSymbolList.this.getTag(symPos);
					}

					public Symbol setSymbol(Symbol sym) {
						return TaggedSymbolList.this.set(symPos, sym);
					}

					public T setTag(T tag) {
						return TaggedSymbolList.this.setTag(symPos, tag);
					}
				};
			}

			public int previousIndex() {
				return this.pos - 1 + offset;
			}

			public void set(TaggedSymbol<T> e) {
				TaggedSymbolList.this.set(this.pos, e.getSymbol(), e.getTag());
			}

			public boolean hasNext() {
				return this.pos < size();
			}

			public TaggedSymbol<T> next() {
				if (!this.removed) {
					this.pos++;
				}
				this.removed = false;
				return new TaggedSymbol<T>() {
					
					private static final long serialVersionUID = 1L;

					private final int symPos = pos;

					public Symbol getSymbol() {
						return TaggedSymbolList.this.get(symPos);
					}

					public T getTag() {
						return TaggedSymbolList.this.getTag(symPos);
					}

					public Symbol setSymbol(Symbol sym) {
						return TaggedSymbolList.this.set(symPos, sym);
					}

					public T setTag(T tag) {
						return TaggedSymbolList.this.setTag(symPos, tag);
					}
				};
			}

			public void remove() {
				if (!this.removed) {
					TaggedSymbolList.this.remove(this.pos);
					this.removed = true;
				}
			}
		};
	}
}
