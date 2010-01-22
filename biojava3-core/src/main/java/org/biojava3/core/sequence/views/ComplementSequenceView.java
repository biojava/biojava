package org.biojava3.core.sequence.views;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.template.AbstractSequenceHoldingSequenceView;
import org.biojava3.core.sequence.template.NucleotideCompoundInterface;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceView;

public class ComplementSequenceView<C extends NucleotideCompoundInterface> extends AbstractSequenceHoldingSequenceView<C> {

  public ComplementSequenceView(Sequence<C> sequence) {
    super(sequence);
  }

  @SuppressWarnings("unchecked")
  public List<C> getAsList() {
    List<C> list = new ArrayList<C>(getLength());
    for(C c: getViewedSequence()) {
      list.add((C)c.getComplement());
    }
    return list;
  }

  @SuppressWarnings("unchecked")
  public C getCompoundAt(int position) {
    return (C)super.getCompoundAt(position).getComplement();
  }

  @SuppressWarnings("unchecked")
  public int getIndexOf(C compound) {
    return super.getIndexOf((C)compound.getComplement());
  }

  @SuppressWarnings("unchecked")
  public int getLastIndexOf(C compound) {
    return super.getLastIndexOf((C)compound.getComplement());
  }

  public String getString() {
    StringBuilder b = new StringBuilder(getLength());
    for(C c: this) {
      b.append(c.getComplement().toString());
    }
    return b.toString();
  }

  public SequenceView<C> getSubSequence(final int start, final int end) {
    return new ComplementSequenceView<C>(getViewedSequence()) {
      public int getStart() {
        return start;
      }
      public int getEnd() {
        return end;
      }
    };
  }

  public Iterator<C> iterator() {
    return new Iterator<C>() {
      private final Iterator<C> iterator = getAsList().iterator();
      public boolean hasNext() {
        return iterator.hasNext();
      }
      @SuppressWarnings("unchecked")
      public C next() {
        return (C)iterator.next().getComplement();
      }
      public void remove() {
        iterator.remove();
      }
    };
  }
}
