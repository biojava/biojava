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

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.loader.ArrayListProxySequenceReader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;

/**
 * A helper class that allows different ways to read a string and create a DNA sequence. Used in FastaReaderHelper
 * and probably a layer that isn't needed
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class DNASequenceCreator implements
    SequenceCreatorInterface<NucleotideCompound> {

  private final CompoundSet<NucleotideCompound> compoundSet;

  /**
   *
   * @param compoundSet
   */
  public DNASequenceCreator(CompoundSet<NucleotideCompound> compoundSet) {
    this.compoundSet = compoundSet;
  }

/**
 *
 * @param sequence The Sequence from a String
 * @param index Currently not used
 * @return
 */
  public AbstractSequence<NucleotideCompound> getSequence(String sequence,
      long index) {
    return new DNASequence(sequence, compoundSet);
  }
/**
 *
 * @param proxyLoader The Sequence from a ProxySequenceReader
 * @param index Currently not used
 * @return
 */
  public AbstractSequence<NucleotideCompound> getSequence(
      ProxySequenceReader<NucleotideCompound> proxyLoader, long index) {
    return new DNASequence(proxyLoader, compoundSet);
  }

  /**
   *
   * @param list
   * @return
   */
  public AbstractSequence<NucleotideCompound> getSequence(
      List<NucleotideCompound> list) {
    ArrayListProxySequenceReader<NucleotideCompound> store = new ArrayListProxySequenceReader<NucleotideCompound>();
    store.setCompoundSet(compoundSet);
    store.setContents(list);
    return new DNASequence(store);
  }
}
