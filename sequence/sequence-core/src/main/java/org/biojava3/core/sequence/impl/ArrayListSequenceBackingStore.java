package org.biojava3.core.sequence.impl;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.AbstractSequenceView;
import org.biojava3.core.sequence.Compound;
import org.biojava3.core.sequence.CompoundNotFoundError;
import org.biojava3.core.sequence.CompoundSet;
import org.biojava3.core.sequence.Sequence;
import org.biojava3.core.sequence.SequenceBackingStore;
import org.biojava3.core.sequence.SequenceView;

public class ArrayListSequenceBackingStore<C extends Compound> implements
		SequenceBackingStore<C> {
	
	private CompoundSet<C> compoundSet;
	private List<C> parsedCompounds = new ArrayList<C>();

	public CharSequence getAsCharSequence() {
		// TODO Optimise/cache.
		StringBuilder builder = new StringBuilder();
		for (C compound : parsedCompounds) {
			builder.append(this.compoundSet.getCharSequenceForCompound(compound));
		}
		return builder;
	}

	public List<C> getAsList() {
		return this.parsedCompounds;
	}

	public C getCompoundAt(int position) {
		return this.parsedCompounds.get(position-1);
	}

	public int getIndexOf(C compound) {
		return this.parsedCompounds.indexOf(compound)+1;
	}

	public int getLastIndexOf(C compound) {
		return this.parsedCompounds.lastIndexOf(compound)+1;
	}

	public int getLength() {
		return this.parsedCompounds.size();
	}

	public Iterator<C> iterator() {
		return this.parsedCompounds.iterator();
	}

	public void setCompoundSet(CompoundSet<C> compoundSet) {
		this.compoundSet = compoundSet;
	}

	public void setContents(CharSequence charSequence) {
		// Horrendously inefficient - pretty much the way the old BJ did things. 
		// TODO Should be optimised.
		this.parsedCompounds.clear();
		for (int i = 0; i < charSequence.length(); ) {
			CharSequence compoundStr = null;
			C compound = null;
			for (int compoundStrLength = 1; compound==null && compoundStrLength<=compoundSet.getMaxSingleCompoundCharSequenceLength(); compoundStrLength++) {
				compoundStr = charSequence.subSequence(i,i+compoundStrLength);
				compound = compoundSet.getCompoundForCharSequence(compoundStr);
			}
			if (compound==null) {
				throw new CompoundNotFoundError(compoundStr);
			} else {
				i+=compoundStr.length();
			}
			this.parsedCompounds.add(compound);
		}
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
				return ArrayListSequenceBackingStore.this;
			}
		};
	}

}
