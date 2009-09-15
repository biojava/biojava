/*
 * RichSequenceHandler.java
 *
 * Created on March 7, 2006, 3:12 PM
 */

package org.biojavax.bio.seq;

import java.util.Iterator;
import java.util.List;

import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;

/**
 * An interface for classes that know how to handle subsequence operations.
 * Implementations may be optimized so that they perform more efficiently in
 * certain conditions. For example a subsequence operation on a huge BioSQL
 * backed <code>RichSequence</code> could be optimized so that the operation
 * is performed more efficiently than dragging the whole sequence to memory and
 * then doing the operation.
 *
 * Implementations of <code>RichSequence</code> should generally delegate
 * <code>symbolAt(int index)</code>, <code>subStr(int start, int end)</code>,
 * <code>subList(int start, int end)</code> and subSequence(int start, int end)
 * to some implementation of this interface.
 *
 * @author Mark Schreiber
 * @author Richard Holland
 * @since 1.5
 */
public interface RichSequenceHandler {
    
  /**
   * Apply an edit to the Sequence as specified by the edit object.
   *
   * <h2>Description</h2>
   *
   * <p>
   * All edits can be broken down into a series of operations that change
   * contiguous blocks of the sequence. This represent a one of those operations.
   * </p>
   *
   * <p>
   * When applied, this Edit will replace 'length' number of symbols starting a
   * position 'pos' by the SymbolList 'replacement'. This allow to do insertions
   * (length=0), deletions (replacement=SymbolList.EMPTY_LIST) and replacements
   * (length>=1 and replacement.length()>=1).
   * </p>
   *
   * <p>
   * The pos and pos+length should always be valid positions on the SymbolList
   * to:
   * <ul>
   * <li>be edited (between 0 and symL.length()+1).</li>
   * <li>To append to a sequence, pos=symL.length()+1, pos=0.</li>
   * <li>To insert something at the beginning of the sequence, set pos=1 and
   * length=0.</li>
   * </ul>
   * </p>
   *
   * <h2>Examples</h2>
   *
   * <code><pre>
   * RichSequence seq = //code to initialize RichSequence
   * System.out.println(seq.seqString());
   *
   * // delete 5 bases from position 4
   * Edit ed = new Edit(4, 5, SymbolList.EMPTY_LIST);
   * seq.edit(ed);
   * System.out.println(seq.seqString());
   *
   * // delete one base from the start
   * ed = new Edit(1, 1, SymbolList.EMPTY_LIST);
   * seq.edit(ed);
   *
   * // delete one base from the end
   * ed = new Edit(seq.length(), 1, SymbolList.EMPTY_LIST);
   * seq.edit(ed);
   * System.out.println(seq.seqString());
   *
   * // overwrite 2 bases from position 3 with "tt"
   * ed = new Edit(3, 2, DNATools.createDNA("tt"));
   * seq.edit(ed);
   * System.out.println(seq.seqString());
   *
   * // add 6 bases to the start
   * ed = new Edit(1, 0, DNATools.createDNA("aattgg");
   * seq.edit(ed);
   * System.out.println(seq.seqString());
   *
   * // add 4 bases to the end
   * ed = new Edit(seq.length() + 1, 0, DNATools.createDNA("tttt"));
   * seq.edit(ed);
   * System.out.println(seq.seqString());
   *
   * // full edit
   * ed = new Edit(3, 2, DNATools.createDNA("aatagaa");
   * seq.edit(ed);
   * System.out.println(seq.seqString());
   * </pre></code>
   *
   * @param edit the Edit to perform
   * @throws IndexOutOfBoundsException if the edit does not lie within the
   *         SymbolList
   * @throws IllegalAlphabetException if the SymbolList to insert has an
   *         incompatible alphabet
   * @throws ChangeVetoException  if either the SymboList does not support the
   *         edit, or if the change was vetoed
   */
    public void edit(RichSequence seq, Edit edit) throws IndexOutOfBoundsException, IllegalAlphabetException, ChangeVetoException;
    
  /**
   * Return the symbol at index, counting from 1.
   *
   * @param index the offset into this SymbolList
   * @return  the Symbol at that index
   * @throws IndexOutOfBoundsException if index is less than 1, or greater than
   *                                   the length of the symbol list
   */
    public Symbol symbolAt(RichSequence seq, int index) throws IndexOutOfBoundsException;
    
   /**
   * Returns a List of symbols.
   * <p>
   * This should be an immutable list of symbols or a copy.
   *
   * @return  a List of Symbols
   */
    public List toList(RichSequence seq);
    
  /**
   * Return a region of this sequence as a String.
   * <p>
   * This should use the same rules as seqString.
   *
   * @param start  the first symbol to include
   * @param end the last symbol to include
   * @return the string representation
   * @throws IndexOutOfBoundsException if either start or end are not within the
   *         SymbolList
   */
    public String subStr(RichSequence seq, int start, int end) throws IndexOutOfBoundsException;
    
  /**
   * Return a new SymbolList for the symbols start to end inclusive.
   * <p>
   * The resulting SymbolList will count from 1 to (end-start + 1) inclusive, and
   * refer to the symbols start to end of the original sequence.
   *
   * @param start the first symbol of the new SymbolList
   * @param end the last symbol (inclusive) of the new SymbolList
   */
    public SymbolList subList(RichSequence seq, int start, int end) throws IndexOutOfBoundsException;
    
   /**
   * Stringify this Sequences.
   * <p>
   * It is expected that this will use the symbol's token to render each
   * symbol. It should be parsable back into a SymbolList using the default
   * token parser for this alphabet.
   *
   * @return  a string representation of the symbol list
   */
    public String seqString(RichSequence seq);
    
   /**
   * An Iterator over all Symbols in this SymbolList.
   * <p>
   * This is an ordered iterator over the Symbols. It cannot be used
   * to edit the underlying symbols.
   *
   * @return  an iterator
   */
    public Iterator iterator(RichSequence seq);
}
