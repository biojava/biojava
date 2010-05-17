package org.biojava3.core.sequence;

import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Bare bones version of the Sequence object to be used sparingly. You should
 * really use a specialised version of {@link Sequence} which describes
 * your domain.
 */
public class BasicSequence<C extends Compound> extends AbstractSequence<C> {

  public BasicSequence(String sequence, CompoundSet<C> compoundSet) {
    super(sequence, compoundSet);
  }

  public BasicSequence(ProxySequenceReader<C> reader) {
    super(reader, reader.getCompoundSet());
  }
}
