package org.biojava3.core.sequence.edits;

import org.biojava3.core.sequence.BasicSequence;
import org.biojava3.core.sequence.storage.JoiningSequenceReader;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

public interface Edit<C extends Compound> {

  Sequence<C> edit(Sequence<C> sequence);

  public static abstract class AbstractEdit<C extends Compound> implements Edit<C> {

    public Sequence<C> edit(Sequence<C> editingSequence) {
      Sequence<C> targetSequence = getSequence(editingSequence);
      return editAction(editingSequence, targetSequence);
    }

    @SuppressWarnings("unchecked")
    protected Sequence<C> editAction(Sequence<C> editingSequence,
        Sequence<C> targetSequence) {
      if(isEditAtStart()) {
        Sequence<C> threePrime = getThreePrime(editingSequence);
        return new JoiningSequenceReader<C>(targetSequence, threePrime);
      }
      else if(isEditAtEnd(editingSequence)) {
        Sequence<C> fivePrime = getFivePrime(editingSequence);
        return new JoiningSequenceReader<C>(fivePrime, targetSequence);
      }
      else {
        Sequence<C> fivePrime = getFivePrime(editingSequence);
        Sequence<C> threePrime = getThreePrime(editingSequence);
        return new JoiningSequenceReader<C>(fivePrime, targetSequence, threePrime);
      }
    }

    private int start;
    private int end = -1;

    private String stringSequence;
    private Sequence<C> sequence;

    public AbstractEdit(int start) {
      this.start = start;
    }

    public AbstractEdit(int start, int end) {
      this.start = start;
      this.end = end;
    }

    protected void setStringSequence(String stringSequence) {
      this.stringSequence = stringSequence;
    }

    protected void setSequence(Sequence<C> sequence) {
      this.sequence = sequence;
    }

    /**
     * Returns the Sequence which is our edit. As part of this code
     * we actually
     *
     * @param editingSequence Asked for in-case we need to do String to
     * Sequence conversion so we need a CompoundSet which is given
     * by the Sequence we are editing
     * @return The Sequence<C> object we wish to insert
     */
    public Sequence<C> getSequence(Sequence<C> editingSequence) {
      if(sequence == null && stringSequence != null) {
        sequence = new BasicSequence<C>(
            stringSequence, editingSequence.getCompoundSet());
        if(end < 0) {
          end = (getStart() + sequence.getLength())-1;
        }
      }
      return sequence;
    }

    public int getStart() {
      return start;
    }

    public int getEnd() {
      return end;
    }

    public int getLength() {
      return (getEnd() - getStart())+1;
    }

    protected boolean isEditAtStart() {
      return (getStart() == 1);
    }

    protected boolean isEditAtEnd(Sequence<C> originalSequence) {
     return (getEnd() == originalSequence.getLength());
    }

    protected void assertEditWithinBounds(Sequence<C> originalSequence) {
      int length = getLength();
      int sequenceLength = originalSequence.getLength();
      if(length < 0) {
        throw new IndexOutOfBoundsException("Length of edit is less than 0; cannot remove a negative length");
      }
      if(getStart() < 1) {
        throw new IndexOutOfBoundsException("Edit start is "+getStart()+" is less than 1");
      }
      if(getEnd() > sequenceLength) {
        throw new IndexOutOfBoundsException("Edit end is "+getEnd()+" is greater than "+sequenceLength);
      }
    }

    protected Sequence<C> getFivePrime(Sequence<C> sequence) {
      int start = getStart()-1;
      if (start < 1) start = 1;
      return sequence.getSubSequence(1, start);
    }

    protected Sequence<C> getThreePrime(Sequence<C> sequence) {
      return sequence.getSubSequence(getEnd()+1, sequence.getLength());
    }
  }

  public static class Delete<C extends Compound> extends AbstractEdit<C> {

    public Delete(int position) {
      this(position, position);
    }

    public Delete(int start, int end) {
      super(start, end);
      setStringSequence("");
    }
  }

  public static class Insert<C extends Compound> extends AbstractEdit<C> {

    public Insert(String sequence, int position) {
      super(position, position);
      setStringSequence(sequence);
    }

    public Insert(Sequence<C> sequence, int position) {
      super(position, position);
      setSequence(sequence);
    }

    /**
     * Unique situation in insertions where the start of the 3' end
     * is not the end of the edit plus 1.
     */
    @Override
    protected Sequence<C> getThreePrime(Sequence<C> sequence) {
      return sequence.getSubSequence(getEnd(), sequence.getLength());
    }

    /**
     * Small override because when we deal with 3 prime edits we still have
     * 3 components to glue together.
     */
    @Override
    @SuppressWarnings("unchecked")
    protected Sequence<C> editAction(Sequence<C> editingSequence,
        Sequence<C> targetSequence) {
      if(isEditAtEnd(editingSequence)) {
        Sequence<C> fivePrime = getFivePrime(editingSequence);
        Sequence<C> threePrime = getThreePrime(editingSequence);
        return new JoiningSequenceReader<C>(fivePrime, targetSequence, threePrime);
      }
      else {
        return super.editAction(editingSequence, targetSequence);
      }
    }
  }

  public static class Substitute<C extends Compound> extends AbstractEdit<C> {

    public Substitute(String sequence, int position) {
      super(position);
      setStringSequence(sequence);
    }

    public Substitute(Sequence<C> sequence, int position) {
      super(position, (position + sequence.getLength())-1 );
      setSequence(sequence);
    }
  }

}
