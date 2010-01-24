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
package org.biojava3.core.sequence.storage;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.template.AbstractSequenceView;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceBackingStore;
import org.biojava3.core.sequence.template.SequenceView;

public class ArrayListSequenceBackingStore<C extends Compound> implements SequenceBackingStore<C> {

    private CompoundSet<C> compoundSet;
    private List<C> parsedCompounds = new ArrayList<C>();

    public String getSequenceAsString() {
        // TODO Optimise/cache.
        StringBuilder builder = new StringBuilder();
        for (C compound : parsedCompounds) {
            builder.append(this.compoundSet.getStringForCompound(compound));
        }
        return builder.toString();
    }

    public List<C> getAsList() {
        return this.parsedCompounds;
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

    public int getLength() {
        return this.parsedCompounds.size();
    }

    public Iterator<C> iterator() {
        return this.parsedCompounds.iterator();
    }

    public void setCompoundSet(CompoundSet<C> compoundSet) {
        this.compoundSet = compoundSet;
    }

    public CompoundSet<C> getCompoundSet() {
      return compoundSet;
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
