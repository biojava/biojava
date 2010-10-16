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

import org.biojava.bio.Annotation;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;


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
   * @author     Thomas Down
   * @author     Matthew Pocock
   */

class SubList extends AbstractSymbolList implements Serializable {
    private SymbolList parent;
    private int start, end;

    private transient EditTranslater editTranslater = null;
    private transient ChangeForwarder annotationForwarder = null;

    public SubList(SymbolList parent, int start, int end) {
      this.parent = parent;
      this.start = start;
      this.end = end;
    }

    public Alphabet getAlphabet() {
      return parent.getAlphabet();
    }

    public int length() {
      return end - start + 1;
    }

    public Symbol symbolAt(int pos) {
      if (pos < 1 || pos > length()) {
        throw new IndexOutOfBoundsException(
            "Symbol index out of bounds " + length() + ":" + pos
            );
      }
      return parent.symbolAt(pos + start - 1);
    }

    public SymbolList subList(int sstart, int send) {
      if (sstart < 1 || send > length()) {
        throw new IndexOutOfBoundsException(
            "Sublist index out of bounds " + length() + ":" + sstart + "," + send
            );
      }

      if (send < sstart) {
        throw new IndexOutOfBoundsException(
            "Requested end must not be lower than start: start=" + sstart + ", end=" + send
            );
      }

      return new SubList(parent, sstart + start - 1, send + start - 1);
    }


    // fixme: doesn't do range checking on edit object
    public void edit(Edit edit) throws IllegalAlphabetException, ChangeVetoException {
      parent.edit(new Edit(
        edit.pos + this.start - 1, edit.length, edit.replacement
      ));
    }

    protected ChangeSupport getChangeSupport(ChangeType changeType) {
      ChangeSupport cs = super.getChangeSupport(changeType);

      if(
        (SymbolList.EDIT.isMatchingType(changeType) || changeType.isMatchingType(SymbolList.EDIT)) &&
        (editTranslater == null)
      ) {
        editTranslater = new EditTranslater(this, cs, start, end);
        parent.addChangeListener(editTranslater, SymbolList.EDIT);
      }

      if(
        ((changeType == null) || (changeType == Annotation.PROPERTY)) &&
        (annotationForwarder == null)
      ) {
        annotationForwarder =
                new ChangeForwarder.Retyper(this, cs, Annotation.PROPERTY);
        parent.addChangeListener(annotationForwarder, Annotation.PROPERTY);
      }

      return cs;
    }


    public String toString() {
      return super.toString() +
             " start: " + start +
             " end: " + end;
    }
}