/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.io;

import java.util.List;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.loader.ArrayListProxySequenceReader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class DNASequenceCreator implements
    SequenceCreatorInterface<NucleotideCompound> {

  private final CompoundSet<NucleotideCompound> compoundSet;

  public DNASequenceCreator(CompoundSet<NucleotideCompound> compoundSet) {
    this.compoundSet = compoundSet;
  }

  public AbstractSequence<NucleotideCompound> getSequence(String sequence,
      long index) {
    return new DNASequence(sequence, compoundSet);
  }

  public AbstractSequence<NucleotideCompound> getSequence(
      ProxySequenceReader<NucleotideCompound> proxyLoader, long index) {
    return new DNASequence(proxyLoader, compoundSet);
  }

  public AbstractSequence<NucleotideCompound> getSequence(
      List<NucleotideCompound> list) {
    ArrayListProxySequenceReader<NucleotideCompound> store = new ArrayListProxySequenceReader<NucleotideCompound>();
    store.setCompoundSet(compoundSet);
    store.setContents(list);
    return new DNASequence(store);
  }
}
