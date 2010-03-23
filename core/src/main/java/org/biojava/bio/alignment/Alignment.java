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

package org.biojava.bio.alignment;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeType;

/**
 * An alignment containing multiple <span class="type">SymbolList</span>s.
 * <p>
 * The alignment can be thought of as a rectangular array of
 * <span class="type">Symbol</span>s. Each
 * row is indexed by a label and each column is indexed by offset (counting from
 * 1).
 * <p>
 * Alternatively, it can be thought of as a <span class="type">SymbolList</span>
 * where each <span class="type">Symbol</span> is
 * a list of <span class="type">Symbol</span>s in that column.
 * <p>
 * To create gapped alignments, use <span class="type">SymbolList</span>s with
 * gaps. The most flexible way to do this will be to leverage
 * <span class="type">GappedSymbolList</span> objects.
 *
 * @author Matthew Pocock
 * @author Nimesh Singh
 * @since 1.1
 */
public interface Alignment extends SymbolList {
  /**
   * Signals that SymbolLists will be added to or removed from an alignment. The
   * ChangeEvent will record Object[] { label, symbolList } in previous if it is
   * being removed, in current if it is being added and in both if the
   * SymbolList for a given name is swapped.
   */
  public static final ChangeType CONTENT = new ChangeType(
    "The sequences contained in this alignment are being changed",
    "org.biojava.bio.symbol.Alignment",
    "CONTENT"
  );

  /**
   * <p>
   * The list of SymbolLists in the alignment.
   * </p>
   *
   * <p>
   * The index in the list is the same as the index in the alignment.
   * Each SymbolList object will only be in the alignment once. However, a
   * single underlying SymbolList may have more than one view within an
   * alignment, each represented by a different GappedSymbolList.
   * </p>
   *
   * @return  the List of all SymbolLists in the alignment
   */
  List<String> getLabels();

  /**
   * Retrieve a symbol by label and column.
   *
   * @param label the SymbolList to retrieve from
   * @param column  the index of the column to retrieve
   * @return  the symbol in the symbol list associated with the label at the given column
   * @throws NoSuchElementException if there is no row for 'label'
   */
  Symbol symbolAt(Object label, int column)
  throws NoSuchElementException;

  /**
   * Retrieve a single row of the alignment by label.
   *
   * @param label the object from which to retrieve the symbol list
   * @return  a SymbolList that contains each token in a row of the alignment
   * @throws NoSuchElementException if there is no row for 'label'
   */
  SymbolList symbolListForLabel(Object label)
  throws NoSuchElementException;

  /**
   * <p>
   * Make a view onto this alignment.
   * </p>
   *
   * <p>
   * If labels is null, then each label will be kept. Otherwise, only those in
   * labels will be kept.
   * If loc is null, then the entire length of the alignment will be kept.
   * If loc is not null, then only the columns within the location will be kept.
   * </p>
   *
   * @param labels the Set of sequences to include by label
   * @param loc the Location to include
   * @return  a sub Alignment
   * @throws  NoSuchElementException if labels contains any item that is not a label
   */
  Alignment subAlignment(Set<String> labels, Location loc)
  throws NoSuchElementException;

  /**
   * Creates an Iterator over the SymbolLists in the alignment. This should be
   * similar to iterating over the labels and then fetching each SymbolList, but
   * the order is not guaranteed to be the same.
   *
   * @return an Iterator
   */
  Iterator<SymbolList> symbolListIterator();
  
  /**
   * Iterator implementation looping over symbol lists in an alignment using
   * the labels. This is intended for Alignment implementors.
   *
   * @author Matthew Pocock
   */
  public static class SymbolListIterator
  implements Iterator<SymbolList> {
    private final Iterator<String> labIt;
    private final Alignment ali;
    
    /*
     * 
     */
    public SymbolListIterator(Alignment ali) {
      this.ali = ali;
      this.labIt = ali.getLabels().iterator();
    }
    
    /*
     * (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    public boolean hasNext() {
      return labIt.hasNext();
    }
    
    /*
     * (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    public SymbolList next() {
      return ali.symbolListForLabel(labIt.next());
    }
    
    /*
     * (non-Javadoc)
     * @see java.util.Iterator#remove()
     */
    public void remove() {
      labIt.remove();
    }
  }
}

