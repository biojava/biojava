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
 * Factory methods for constructing useful views onto SymbolLists.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.1
 */

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;

/**
 * Tools class for constructing views of <code>SymbolList</code> objects.
 *
 * @since 1.2
 */

public final class SymbolListViews {
  private SymbolListViews() {}

    /**
     * An n-th order view of another SymbolList.
     * <p>
     * In practice, what this means is that you can view a DNA sequence into an
     * overlapping dinucleotide sequence without having to do any work yourself.
     * </p>
     *
     * @param source The underlying SymbolList to view
     * @param order The window size
     */

    public static SymbolList orderNSymbolList(SymbolList source, int order)
        throws IllegalAlphabetException
    {
	if (order == 1)
	    return source;

	return new OrderNSymbolList(source, order);
    }

    /**
     * A view of windows onto another SymbolList.
     * <p>
     * In practice, what this means is that you can view a DNA sequence as codons which
     * do not overlap.
     * </p>
     *
     * @param source The underlying SymbolList to view
     * @param wsize The window size.
     * @throws IllegalArgumentException if the symbollist length isn't an integer multiple of wsize.
     */

    public static SymbolList windowedSymbolList(SymbolList source, int wsize)
        throws IllegalArgumentException
    {
	return new WindowedSymbolList(source, wsize);
    }

    /**
     * A reversed view onto a SymbolList.
     *
     * @param symbols the SymbolList to reverse.
     */

    public static SymbolList reverse(SymbolList symbols) {
	return new ReverseSymbolList(symbols);
    }

    /**
     * Provides a 'translated' view of an underlying SymbolList.
     * <p>
     * This method allows you to translate from one alphabet into another, so
     * for example, you could translate from DNA-triplets into amino-acids. You
     * could also translate from DNA-dinucleotide into the 'twist' structural
     * metric, or any other translation that takes your fancy.
     * </p>
     * <p>
     * The actual mapping from source to view Symbol is encapsulated in a
     * TranslationTable object.
     * </p>
     * <p>
     * The translated SymbolList will be the same length as the source, and each
     * Symbol in the view will correspond to a single Symbol in the source.
     * </p>
     *
     * @param symbols a SymbolList to translate.
     * @param table a translation table for mapping symbols.
     */

    public static SymbolList translate(SymbolList symbols,
				TranslationTable table)
        throws IllegalAlphabetException
    {
	return new TranslatedSymbolList(symbols, table);
    }

    /**
     * Construct an alignment of the SymbolLists contained in the values collection
     * of <code>labelToSymList</code>.
     *
     * @param labelToSymList A Map containing label -> SymbolList mappings
     */

    public static Alignment alignment(Map labelToSymList)
    throws IllegalArgumentException {
      return new SimpleAlignment(labelToSymList);
    }

    /**
     * View a SymbolList over a cross-product Alphabet as an Alignment.
     *
     * @param labels a List of labels, which should be the same length
     *               as the order <code>symList</code>'s Alphabet.
     * @param symList a SymbolList over a cross-product alphabet.
     */


    public static Alignment alignment(List labels, SymbolList symList)
    throws IllegalArgumentException {
      return new SymListAsAlignment(labels, symList);
    }

    /**
     * View a portion of a SymbolList.  Unlike SymbolList.subList, this
     * method is guarenteed to return a view, which will change when
     * the underlying SymbolList is modified.
     *
     * @param parent the SymbolList to view
     * @param start the first index to include in the view
     * @param end the last index to include in the view
     * @throws IllegalArgumentException if the start or end points fall outside the parent SymbolList.
     * @since 1.4
     */

    public static SymbolList subList(SymbolList parent, int start, int end)
        throws IllegalArgumentException
    {
        if (start < 1 || end > parent.length()) {
            throw new IndexOutOfBoundsException(
                "Sublist index out of bounds " + parent.length() + ":" + start + "," + end
            );
        }

        if (end < start) {
            throw new IllegalArgumentException(
                "end must not be lower than start: start=" + start + ", end=" + end
            );
        }
        return new SubList(parent, start, end);
    }

  /**
   * Get a new immutable, empty symbol list with the given alphabet.
   *
   * @since 1.4
   * @param alpha   the Alphabet this symbol list is over
   * @return  a new empty SymbolList
   */
  public static SymbolList emptyList(Alphabet alpha)
  {
    return new EmptySymbolList(alpha);
  }

    private static class SymListAsAlignment
    extends Unchangeable
    implements Alignment {
      private final SymbolList symList;
      private final List labels;

      public SymListAsAlignment(List labels, SymbolList symList) {
        if(labels.size() != symList.getAlphabet().getAlphabets().size()) {
          throw new IllegalArgumentException("There must be one label per symbol list");
        }

        this.labels = Collections.unmodifiableList(new ArrayList(labels));
        this.symList = symList;
      }

      public List getLabels() {
        return labels;
      }

      public SequenceIterator sequenceIterator() {
        throw new UnsupportedOperationException("This method sucks");
      }

      public Iterator symbolListIterator() {
        return new Alignment.SymbolListIterator(this);
      }

      public Symbol symbolAt(Object label, int column) {
        BasisSymbol sym = (BasisSymbol) symList.symbolAt(column);
        return (Symbol) sym.getSymbols().get(labels.indexOf(label));
      }

      public SymbolList symbolListForLabel(Object label) {
        return new IndexedSymbolList(symList, labels.indexOf(label));
      }

      public Alphabet getAlphabet() {
        return symList.getAlphabet();
      }

      public Iterator iterator() {
        return symList.iterator();
      }

      public int length() {
        return symList.length();
      }

      public String seqString() {
        return symList.seqString();
      }

      public SymbolList subList(int start, int end)
      throws IndexOutOfBoundsException {
        return symList.subList(start, end);
      }

      public String subStr(int start, int end)
      throws IndexOutOfBoundsException {
        return symList.subStr(start, end);
      }

      public void edit(Edit edit)
      throws
        IndexOutOfBoundsException,
        IllegalAlphabetException,
        ChangeVetoException
      {
        symList.edit(edit);
      }

      public List toList() {
        return symList.toList();
      }

      public Symbol symbolAt(int indx)
      throws IndexOutOfBoundsException {
        return symList.symbolAt(indx);
      }

      public Alignment subAlignment(Set labels, Location loc) {
        throw new UnsupportedOperationException("Fixme: this needs to be implemented");
      }
    }

    private static class IndexedSymbolList
    extends AbstractSymbolList {
      private final int indx;
      private final SymbolList symList;

      public IndexedSymbolList(SymbolList symList, int indx)
      throws IllegalArgumentException {
        if(indx >= symList.getAlphabet().getAlphabets().size()) {
          throw new IllegalArgumentException("index too high");
        }

        this.indx = indx;
        this.symList = symList;
      }

      public Alphabet getAlphabet() {
        return (Alphabet) symList.getAlphabet().getAlphabets().get(indx);
      }

      public int length() {
        return symList.length();
      }

      public Symbol symbolAt(int indx)
      throws IndexOutOfBoundsException {
        return (Symbol) ((BasisSymbol) symList.symbolAt(indx)).getSymbols().get(this.indx);
      }
    }
}
