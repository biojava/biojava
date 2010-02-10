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

package org.biojava.bio.symbol;

/**
 * This extends SymbolList with API for manipulating, inserting and deleting
 * gaps.
 * <p>
 * You could make a SymbolList that contains gaps directly. However, this leaves
 * you with a nasty problem if you wish to support gap-edit operations. Also,
 * the original SymbolList must either be coppied or lost.
 * <p>
 * GappedSymbolList solves these problems. It will maintain a data-structure
 * that places gaps. You can add and remove the gaps by using the public API.
 * <p>
 * For gap-insert operations, the insert index is the position that will become
 * a gap. The symbol currently there will move to a higher index. To insert
 * leading gaps, add gaps at index 1. To insert trailing gaps, add gaps at index
 * length+1.
 * 
 * @author Matthew Pocock
 * @since 1.3
 */
public interface GappedSymbolList extends SymbolList {
	/**
	 * Return the underlying (ungapped) SymbolList.
	 * 
	 * @since 1.4
	 */
	public SymbolList getSourceSymbolList();

	/**
	 * Coordinate conversion from view to source.
	 * <p>
	 * If the index can be projected onto the source, the index it projects onto
	 * is returned. If it falls within a gap, then the index of the first symbol
	 * after the run of gaps is negated and returned. If the index is after the
	 * last block of symbols (and therefore in the trailing list of gaps), then
	 * it returns -(length + 1).
	 * 
	 * @param indx
	 *            the index to project
	 * @return the position of indx projected from view to source
	 * @throws IndexOutOfBoundsException
	 *             if indx is not a valid view index
	 */
	public int viewToSource(int indx) throws IndexOutOfBoundsException;

	/**
	 * Coordinate conversion from source to view.
	 * 
	 * @param indx
	 *            the index to project
	 * @return the position of indx projected from source to view
	 * @throws IndexOutOfBoundsException
	 *             if indx is not a valid source index
	 */
	public int sourceToView(int indx) throws IndexOutOfBoundsException;

	/**
	 * Add a single gap at pos within the view coordintates.
	 * <p>
	 * this.symbolAt(pos) will then return gap. Adding a gap at 1 will prepend
	 * gaps. Adding a gap at (length+1) will append a gap.
	 * 
	 * @param pos
	 *            the position to add a gap before
	 * @throws IndexOutOfBoundsException
	 *             if pos is not within 1->length+1
	 */
	public void addGapInView(int pos) throws IndexOutOfBoundsException;

	/**
	 * Add length gaps at pos within the view coordinates.
	 * <p>
	 * this.symbolAt(i) will then return gap for i = (pos .. pos+count-1).
	 * Adding gaps at 1 will prepend gaps. Adding gaps at (length+1) will append
	 * gaps.
	 * 
	 * @param pos
	 *            the position to add a gap before
	 * @param length
	 *            the number of gaps to insert
	 * @throws IndexOutOfBoundsException
	 *             if pos is not within 1->length+1
	 */
	public void addGapsInView(int pos, int length)
			throws IndexOutOfBoundsException;

	/**
	 * Add a gap at pos within the source coordinates.
	 * 
	 * @param pos
	 *            where to add the gap
	 * @throws IndexOutOfBoundsException
	 *             if pos is not within 1->source.length()
	 */
	public void addGapInSource(int pos) throws IndexOutOfBoundsException;

	/**
	 * Add length gaps at pos within the source coordinates.
	 * 
	 * @param pos
	 *            where to add the gap
	 * @param length
	 *            how many gaps to add
	 * @throws IndexOutOfBoundsException
	 *             if pos is not within 1->source.length()
	 */
	public void addGapsInSource(int pos, int length)
			throws IndexOutOfBoundsException;

	/**
	 * Remove a single gap at position pos in this GappedSymbolList.
	 * 
	 * @param pos
	 *            where to remove the gap
	 * @throws IndexOutOfBoundsException
	 *             if pos is not within 1..length
	 * @throws IllegalSymbolException
	 *             if the symbol at pos is not a gap
	 */
	public void removeGap(int pos) throws IndexOutOfBoundsException,
			IllegalSymbolException;

	/**
	 * Remove some gaps at position pos in this GappedSymbolList.
	 * 
	 * @param pos
	 *            where to remove the gaps
	 * @param length
	 *            how many to remove
	 * @throws IndexOutOfBoundsException
	 *             if pos is not within 1..length() or if some of the Symbols
	 *             within pos->(pos+length-1) are not gap Symbols
	 * @throws IllegalSymbolException
	 *             if the symbol at pos is not a gap
	 */
	public void removeGaps(int pos, int length)
			throws IndexOutOfBoundsException, IllegalSymbolException;

	/**
	 * Return the index of the first Symbol that is not a Gap character.
	 * <p>
	 * All symbols before firstNonGap are leading gaps. firstNonGap is
	 * effectively the index in the view of symbol 1 in the original sequence.
	 * 
	 * @return the index of the first character not to be a gap
	 */
	public int firstNonGap();

	/**
	 * Return the index of the last Symbol that is not a Gap character.
	 * <p>
	 * All symbols after lastNonGap untill length are trailing gaps. lastNonGap
	 * is effectively the index in the view of symbol length in the original
	 * sequence.
	 * 
	 * @return the index of the last character not to be a gap
	 */
	public int lastNonGap();

	/**
	 * Get a Location that contains exactly those positions that are not gaps.
	 * 
	 * <p>
	 * This will be a Location that contains every symbol in the underlying
	 * ungapped sequence. Every symbol not in the Location is not from the
	 * underlying sequence and is a gap.
	 * </p>
	 * 
	 * @return a new Location that contains all non-gap symbols
	 */
	public Location getUngappedLocation();

}
