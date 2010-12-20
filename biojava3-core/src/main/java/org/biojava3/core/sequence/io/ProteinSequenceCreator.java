/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 01-21-2010
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
 * Used to create a ProteinSequence from a String to allow for details
 * about the location of the sequence etc.
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class ProteinSequenceCreator implements
    SequenceCreatorInterface<AminoAcidCompound> {

  private CompoundSet<AminoAcidCompound> compoundSet;
/**
 *
 * @param compoundSet
 */
  public ProteinSequenceCreator(CompoundSet<AminoAcidCompound> compoundSet) {
    this.compoundSet = compoundSet;
  }
/**
 *
 * @param sequence
 * @param index not used in this implementation
 * @return
 */
  public AbstractSequence<AminoAcidCompound> getSequence(String sequence,
      long index) {
    return new ProteinSequence(sequence, compoundSet);
  }
/**
 *
 * @param list
 * @return
 */
  public AbstractSequence<AminoAcidCompound> getSequence(
      List<AminoAcidCompound> list) {
    ArrayListProxySequenceReader<AminoAcidCompound> store = new ArrayListProxySequenceReader<AminoAcidCompound>();
    store.setCompoundSet(compoundSet);
    store.setContents(list);
    return new ProteinSequence(store);
  }
/**
 *
 * @param proxyLoader
 * @param index not used in this implementation
 * @return
 */
  public AbstractSequence<AminoAcidCompound> getSequence(
      ProxySequenceReader<AminoAcidCompound> proxyLoader, long index) {
    return new ProteinSequence(proxyLoader, compoundSet);
  }
}
