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
