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
	
	public String getString() {
		// TODO Optimise.
		return this.getViewedSequence().getString().substring(this.getStart()-1, this.getEnd());
	}

	public List<C> getAsList() {
		// TODO Optimise.
		return this.getViewedSequence().getAsList().subList(this.getStart()-1, this.getEnd());
	}

	public C getCompoundAt(int position) {
		return this.getViewedSequence().getCompoundAt(this.getStart()+position);
	}

	public int getIndexOf(C compound) {
		return this.getViewedSequence().getIndexOf(compound)+this.getStart();
	}

	public int getLastIndexOf(C compound) {
		return this.getViewedSequence().getLastIndexOf(compound)+this.getStart();
	}

	public int getLength() {
		return (this.getEnd()-this.getStart()) + 1;
	}

	public SequenceView<C> getSubSequence(final int start, final int end) {
		return new AbstractSequenceView<C>() {

			public int getEnd() {
				return end;
			}

			public int getStart() {
				return start;
			}

			public Sequence<C> getViewedSequence() {
				return AbstractSequenceView.this;
			}
		};
	}

	public Iterator<C> iterator() {
		// TODO Optimise.
		return this.getAsList().iterator();
	}
}
