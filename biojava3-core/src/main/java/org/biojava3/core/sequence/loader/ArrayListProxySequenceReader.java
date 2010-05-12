package org.biojava3.core.sequence.loader;

import org.biojava3.core.sequence.storage.ArrayListSequenceReader;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.ProxySequenceReader;

public class ArrayListProxySequenceReader<C extends Compound>
  extends ArrayListSequenceReader<C> implements ProxySequenceReader<C>{

}
