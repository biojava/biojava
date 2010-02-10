package org.biojava3.core.sequence.loader;

import org.biojava3.core.sequence.storage.ArrayListSequenceBackingStore;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.SequenceProxyLoader;

public class SequenceArrayListProxyLoader<C extends Compound>
  extends ArrayListSequenceBackingStore<C> implements SequenceProxyLoader<C>{

}
