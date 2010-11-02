package org.biojava3.core.sequence.views;

import org.biojava3.core.sequence.template.ComplementCompound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceMixin;
import org.biojava3.core.sequence.template.SequenceProxyView;

/**
 * For a given sequence this class will create a view over the top of it
 * and for every request the code will return the complement of the underlying
 * base e.g. base A will become base T
 *
 * @author Andy Yates
 * @param <C> Must be a subtype of @{link ComplementCompound} since
 * only those support complements
 */
public class ComplementSequenceView<C extends ComplementCompound> extends SequenceProxyView<C> {

    public ComplementSequenceView(Sequence<C> sequence) {
        super(sequence);
    }

    @Override
    public String getSequenceAsString() {
        return SequenceMixin.toString(this);
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
}
