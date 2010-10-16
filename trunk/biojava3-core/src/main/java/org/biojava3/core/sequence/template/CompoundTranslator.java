package org.biojava3.core.sequence.template;

import java.util.List;

public interface CompoundTranslator<F extends Compound, T extends Compound> {

  T translate(F fromCompound);

  List<T> translateMany(F fromCompound);

  Sequence<T> createSequence(Sequence<F> originalSequence);

  List<Sequence<T>> createSequences(Sequence<F> originalSequence);

}
