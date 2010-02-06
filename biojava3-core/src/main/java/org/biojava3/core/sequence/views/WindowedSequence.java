package org.biojava3.core.sequence.views;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

public class WindowedSequence<C extends Compound> implements Iterable<List<C>>{

  private final Sequence<C> sequence;
  private final int windowSize;

  public WindowedSequence(Sequence<C> sequence, int windowSize) {
    this.sequence = sequence;
    this.windowSize = windowSize;
  }

  public int getWindowSize() {
    return windowSize;
  }

  public Sequence<C> getBackingSequence() {
    return sequence;
  }

  /**
   * Calculates start index according to the equation start = ( (index-1) -
   * windowSize) +1
   */
  protected int toStartIndex(int index) {
    return ((index - 1) * getWindowSize()) + 1;
  }

  public int getLength() {
    return getBackingSequence().getLength() / getWindowSize();
  }

  public List<C> getCompounds(int index) {
    int start = toStartIndex(index);
    int window = getWindowSize();
    List<C> output = new ArrayList<C>(3);
    for(int i = 0; i < window; i++) {
      output.add(getBackingSequence().getCompoundAt(start+i));
    }
    return output;
  }

  @Override
  public Iterator<List<C>> iterator() {
    final int start = 1;
    final int end = getLength();
    return new Iterator<List<C>>() {
      private int currentIndex = start;
      @Override
      public boolean hasNext() {
        return currentIndex <= end;
      }
      @Override
      public List<C> next() {
        return WindowedSequence.this.getCompounds(currentIndex++);
      }
      public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a Windowed view");
      }
    };
  }
}
