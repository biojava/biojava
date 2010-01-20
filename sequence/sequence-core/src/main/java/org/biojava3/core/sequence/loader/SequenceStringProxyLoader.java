/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.loader;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.biojava3.core.sequence.storage.ArrayListSequenceBackingStore;
import org.biojava3.core.sequence.template.AbstractSequenceView;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundNotFoundError;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceProxyLoader;
import org.biojava3.core.sequence.template.SequenceView;

/**
 *
 * @author Scooter
 */
public class SequenceStringProxyLoader<C extends Compound> implements SequenceProxyLoader<C> {

    private String sequence;
    private CompoundSet<C> compoundSet;
    private List<C> parsedCompounds = new ArrayList<C>();

    public SequenceStringProxyLoader(String sequence, CompoundSet<C> compoundSet) {
        this.sequence = sequence;
        setCompoundSet(compoundSet);
    }

    public void setCompoundSet(CompoundSet<C> compoundSet) {
        this.compoundSet = compoundSet;
    }

    public void setContents(String sequence) {
        // Horrendously inefficient - pretty much the way the old BJ did things.
        // TODO Should be optimised.
        this.parsedCompounds.clear();
        for (int i = 0; i < sequence.length();) {
            String compoundStr = null;
            C compound = null;
            for (int compoundStrLength = 1; compound == null && compoundStrLength <= compoundSet.getMaxSingleCompoundStringLength(); compoundStrLength++) {
                compoundStr = sequence.substring(i, i + compoundStrLength);
                compound = compoundSet.getCompoundForString(compoundStr);
            }
            if (compound == null) {
                throw new CompoundNotFoundError(compoundStr);
            } else {
                i += compoundStr.length();
            }
            this.parsedCompounds.add(compound);
        }
    }

    public int getLength() {
        return this.parsedCompounds.size();
    }

    public C getCompoundAt(int position) {
        return this.parsedCompounds.get(position - 1);
    }

    public int getIndexOf(C compound) {
        return this.parsedCompounds.indexOf(compound) + 1;
    }

    public int getLastIndexOf(C compound) {
        return this.parsedCompounds.lastIndexOf(compound) + 1;
    }

    @Override
    public String toString() {
        return getString();
    }

    public String getString() {
        return sequence;
    }

    public List<C> getAsList() {
        return this.parsedCompounds;
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
                return SequenceStringProxyLoader.this;
            }
        };
    }

    public Iterator<C> iterator() {
        return this.parsedCompounds.iterator();
    }
}
