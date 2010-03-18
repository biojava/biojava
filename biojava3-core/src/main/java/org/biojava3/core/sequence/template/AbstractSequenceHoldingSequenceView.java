package org.biojava3.core.sequence.template;

import org.biojava3.core.sequence.AccessionID;



public class AbstractSequenceHoldingSequenceView<C extends Compound> extends AbstractSequenceView<C> {

  private final Sequence<C> sequence;

  public AbstractSequenceHoldingSequenceView(Sequence<C> sequence) {
    this.sequence = sequence;
  }

  public int getEnd() {
    return getViewedSequence().getLength();
  }

  public int getStart() {
    return 1;
  }

  public Sequence<C> getViewedSequence() {
    return sequence;
  }

    @Override
    public AccessionID getAccession() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
