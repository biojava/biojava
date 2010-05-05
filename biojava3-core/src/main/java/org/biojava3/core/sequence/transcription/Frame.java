package org.biojava3.core.sequence.transcription;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.views.ComplementSequenceView;
import org.biojava3.core.sequence.views.ReversedSequenceView;

/**
 * Indicates a way of translating a sequence.
 *
 * @author ayates
 */
public enum Frame {
  ONE(1, false),
  TWO(2, false),
  THREE(3, false),
  REVERSED_ONE(1, true),
  REVERSED_TWO(2, true),
  REVERSED_THREE(3, true);

  private final boolean reverse;
  private final int start;

  private Frame(int start, boolean reverse) {
    this.start = start;
    this.reverse = reverse;
  }

  public static Frame getDefaultFrame() {
    return ONE;
  }

  /**
   * Returns all frames in the forward orientation
   */
  public static Frame[] getForwardFrames() {
    return new Frame[]{ONE,TWO,THREE};
  }

  /**
   * Returns all frames which are in the reverse orientation
   */
  public static Frame[] getReverseFrames() {
    return new Frame[]{REVERSED_ONE,REVERSED_TWO,REVERSED_THREE};
  }

  /**
   * Delegates to {@link Frame#values()}
   */
  public static Frame[] getAllFrames() {
    return Frame.values();
  }

  /**
   * Optionally wraps a Sequence in a reverse complementing view (if the
   * frame is on the reverse strand) and creates a sub sequence view if
   * it is required.
   *
   * If you pass in {@link #ONE} in you will get the same {@link Sequence}
   * returned.
   */
  public <C extends NucleotideCompound> Sequence<C> wrap(Sequence<C> incoming) {
    Sequence<C> reversed;
    if(reverse) {
      reversed = new ComplementSequenceView<C>(new ReversedSequenceView<C>(incoming));
    }
    else {
      reversed = incoming;
    }

    Sequence<C> finalSeq;

    if(start == 1) {
      finalSeq = reversed;
    }
    else {
      finalSeq = reversed.getSubSequence(start, reversed.getLength());
    }

    return finalSeq;
  }

}