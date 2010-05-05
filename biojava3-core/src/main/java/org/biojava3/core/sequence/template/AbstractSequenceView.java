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
 * Created on 01-21-2010
 *
 * @author Richard Holland
 *
 *
 */
package org.biojava3.core.sequence.template;

import java.util.Iterator;
import java.util.List;

public abstract class AbstractSequenceView<C extends Compound> implements SequenceView<C> {

  public List<C> getAsList() {
    return SequenceMixin.toList(this);
  }

  public String getSequenceAsString() {
    return SequenceMixin.toStringBuilder(this).toString();
  }

	public C getCompoundAt(int position) {
		return getViewedSequence().getCompoundAt((getStart()+position)-1);
	}

	public int getIndexOf(C compound) {
		return (getViewedSequence().getIndexOf(compound)+getStart())-1;
	}

	public int getLastIndexOf(C compound) {
		return (getViewedSequence().getLastIndexOf(compound)+getStart())-1;
	}

	public int getLength() {
		return (getEnd()-getStart()) + 1;
	}

	public CompoundSet<C> getCompoundSet() {
	  return getViewedSequence().getCompoundSet();
	}

	public SequenceView<C> getSubSequence(final int start, final int end) {
		return SequenceMixin.createSubSequence(this, start, end);
	}

	public Iterator<C> iterator() {
		return SequenceMixin.createIterator(this);
	}

	public int countCompounds(C... compounds) {
    return SequenceMixin.countCompounds(this, compounds);
  }
}
