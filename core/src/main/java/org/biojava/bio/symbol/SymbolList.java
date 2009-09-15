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

import java.util.Iterator;
import java.util.List;

import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * A sequence of symbols that belong to an alphabet.
 * <p>
 * This uses biological coordinates (1 to length).
 *
 * @author Matthew Pocock
 * @author Mark Schreiber
 * @author Francois Pepin
 */
public interface SymbolList extends Changeable {
  /**
   * Signals that the SymbolList is being edited. The getChange field of the
   * event should contain the SymbolList.Edit object describing the change.
   */
  public static final ChangeType EDIT = new ChangeType(
    "the SymbolList has been edited",
    "org.biojava.bio.symbol.SymbolList",
    "EDIT"
  );
  
  /**
   * The alphabet that this SymbolList is over.
   * <p>
   * Every symbol within this SymbolList is a member of this alphabet.
   * <code>alphabet.contains(symbol) == true</code>
   * for each symbol that is within this sequence.
   *
   * @return  the alphabet
   */
  Alphabet getAlphabet();
  
  /**
   * The number of symbols in this SymbolList.
   *
   * @return  the length
   */
  int length();

  /**
   * Return the symbol at index, counting from 1.
   *
   * @param index the offset into this SymbolList
   * @return  the Symbol at that index
   * @throws IndexOutOfBoundsException if index is less than 1, or greater than
   *                                   the length of the symbol list
   */
  Symbol symbolAt(int index) throws IndexOutOfBoundsException;
  
  /**
   * Returns a List of symbols.
   * <p>
   * This is an immutable list of symbols. Do not edit it.
   *
   * @return  a List of Symbols
   */
  List toList();
  
  /**
   * An Iterator over all Symbols in this SymbolList.
   * <p>
   * This is an ordered iterator over the Symbols. It cannot be used
   * to edit the underlying symbols.
   *
   * @return  an iterator
   */
  Iterator iterator();
  
  /**
   * Return a new SymbolList for the symbols start to end inclusive.
   * <p>
   * The resulting SymbolList will count from 1 to (end-start + 1) inclusive, and
   * refer to the symbols start to end of the original sequence.
   *
   * @param start the first symbol of the new SymbolList
   * @param end the last symbol (inclusive) of the new SymbolList
   */
  SymbolList subList(int start, int end) throws IndexOutOfBoundsException;
    
  /**
   * Stringify this symbol list.
   * <p>
   * It is expected that this will use the symbol's token to render each
   * symbol. It should be parsable back into a SymbolList using the default
   * token parser for this alphabet.
   *
   * @return  a string representation of the symbol list
   */
  String seqString();
  
  /**
   * Return a region of this symbol list as a String.
   * <p>
   * This should use the same rules as seqString.
   *
   * @param start  the first symbol to include
   * @param end the last symbol to include
   * @return the string representation
   * @throws IndexOutOfBoundsException if either start or end are not within the
   *         SymbolList
   */
  String subStr(int start, int end) throws IndexOutOfBoundsException;
  
  /**
   * Apply an edit to the SymbolList as specified by the edit object.
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
   * SymbolList seq = DNATools.createDNA("atcaaaaacgctagc");
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
  void edit(Edit edit)
  throws IndexOutOfBoundsException, IllegalAlphabetException,
  ChangeVetoException;
  
  /**
   * A useful object that represents an empty symbol list, to avoid returning
   * null.
   *
   */
  static final SymbolList EMPTY_LIST = new EmptySymbolList();

}
