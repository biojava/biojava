/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.io;

import java.util.List;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.loader.ArrayListProxySequenceReader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class ProteinSequenceCreator implements
    SequenceCreatorInterface<AminoAcidCompound> {

  private CompoundSet<AminoAcidCompound> compoundSet;

  public ProteinSequenceCreator(CompoundSet<AminoAcidCompound> compoundSet) {
    this.compoundSet = compoundSet;
  }

  public AbstractSequence<AminoAcidCompound> getSequence(String sequence,
      long index) {
    return new ProteinSequence(sequence, compoundSet);
  }

  public AbstractSequence<AminoAcidCompound> getSequence(
      List<AminoAcidCompound> list) {
    ArrayListProxySequenceReader<AminoAcidCompound> store = new ArrayListProxySequenceReader<AminoAcidCompound>();
    store.setCompoundSet(compoundSet);
    store.setContents(list);
    return new ProteinSequence(store);
  }

  public AbstractSequence<AminoAcidCompound> getSequence(
      ProxySequenceReader<AminoAcidCompound> proxyLoader, long index) {
    return new ProteinSequence(proxyLoader, compoundSet);
  }
}
