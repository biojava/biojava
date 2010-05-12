package org.biojava3.core.sequence.views;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceProxyView;
import org.biojava3.core.sequence.template.SequenceView;

/**
 * For a given sequence this class will create a view over the top of it
 * and for every request the code will return the complement of the underlying
 * base e.g. base A will become base T
 *
 * @author Andy Yates
 * @param <C> Must be a subtype of @{link ComplementCompound} since
 * only those support complements
 */

//took out interface to get it to compile needs to be fixed
public class ComplementSequenceView<C extends NucleotideCompound> extends SequenceProxyView<C> {

    public ComplementSequenceView(Sequence<C> sequence) {
        super(sequence);
    }

    @SuppressWarnings("unchecked")
    @Override
    public List<C> getAsList() {
        List<C> list = new ArrayList<C>(getLength());
        for (C c : getViewedSequence()) {
            list.add((C) c.getComplement());
        }
        return list;
    }

    @SuppressWarnings("unchecked")
    @Override
    public C getCompoundAt(int position) {
        return (C) super.getCompoundAt(position).getComplement();
    }

    @SuppressWarnings("unchecked")
    @Override
    public int getIndexOf(C compound) {
        return super.getIndexOf((C) compound.getComplement());
    }

    @SuppressWarnings("unchecked")
    @Override
    public int getLastIndexOf(C compound) {
        return super.getLastIndexOf((C) compound.getComplement());
    }

    @Override
    public String getSequenceAsString() {
        StringBuilder b = new StringBuilder(getLength());
        for (C c : this) {
            b.append(c.getComplement().toString());
        }
        return b.toString();
    }

    @Override
    public SequenceView<C> getSubSequence(final Integer start, final Integer end) {
        return new ComplementSequenceView<C>(getViewedSequence()) {

            public Integer getStart() {
                return start;
            }

            public Integer getEnd() {
                return end;
            }
        };
    }

    @Override
    public Iterator<C> iterator() {
        return new Iterator<C>() {

            private final Iterator<C> iterator = getAsList().iterator();

            public boolean hasNext() {
                return iterator.hasNext();
            }

            @SuppressWarnings("unchecked")
            public C next() {
                return (C) iterator.next().getComplement();
            }

            public void remove() {
                iterator.remove();
            }
        };
    }
}
