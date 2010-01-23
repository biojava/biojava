package org.biojava3.core.sequence.views;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.template.AbstractSequenceHoldingSequenceView;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceView;

/**
 * For a given sequence this class will return the base at the reversed
 * position i.e. in a sequence of size 10, if you request base 2 you will get
 * back the base at position 9. Sub-views can be made of this class which
 * also respect the reversed calls.
 *
 * @author Andy Yates
 * @param <C> Must be a subtype of @{link Compound}
 */
public class ReversedSequenceView<C extends Compound> extends AbstractSequenceHoldingSequenceView<C> {

  private final int sequenceSize;

  public ReversedSequenceView(Sequence<C> sequence) {
    super(sequence);
    this.sequenceSize = sequence.getLength();
  }

  protected int toIndex(int index) {
    return (sequenceSize-index)+1;
  }

  @Override
  public List<C> getAsList() {
    //TODO better version than reversing please maybe ...
    List<C> l = super.getAsList();
    Collections.reverse(l);
    return l;
  }

  @Override
  public String getString() {
    //TODO This is because of the way the super class does the create; do a better version please!
    StringBuilder b = new StringBuilder(super.getString());
    b.reverse();
    return b.toString();
  }

  public C getCompoundAt(int position) {
    return super.getCompoundAt(toIndex(position));
  }

  public int getIndexOf(C compound) {
    return toIndex(super.getIndexOf(compound));
  }

  public int getLastIndexOf(C compound) {
    return toIndex(super.getLastIndexOf(compound));
  }

  public SequenceView<C> getSubSequence(final int start, final int end) {
    return new ReversedSequenceView<C>(getViewedSequence()) {
      public int getEnd() {
          return toIndex(end);
      }
      public int getStart() {
          return toIndex(start);
      }
    };
  }

  public Iterator<C> iterator() {
    return new Iterator<C>() {
      private int currentIndex = getStart();
      public boolean hasNext() {
        return currentIndex <= getEnd();
      }
      public C next() {
        return getCompoundAt(currentIndex++);
      }
      public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a Sequence from an Iterator");
      }
    };
  }
}
