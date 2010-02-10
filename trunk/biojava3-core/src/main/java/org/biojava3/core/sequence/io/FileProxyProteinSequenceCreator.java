/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import java.io.File;
import java.util.List;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.loader.SequenceFileProxyLoader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.SequenceProxyLoader;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FileProxyProteinSequenceCreator implements
    SequenceCreatorInterface<AminoAcidCompound> {

  CompoundSet<AminoAcidCompound> compoundSet = null;
  File                           fastaFile   = null;

  public FileProxyProteinSequenceCreator(File fastaFile,
      CompoundSet<AminoAcidCompound> compoundSet) {
    this.compoundSet = compoundSet;
    this.fastaFile = fastaFile;
  }

  public AbstractSequence<AminoAcidCompound> getSequence(String sequence,
      long index) {
    SequenceFileProxyLoader<AminoAcidCompound> sequenceFileProxyLoader = new SequenceFileProxyLoader<AminoAcidCompound>(
        fastaFile, new FastaSequenceParser(), index, sequence.length(),
        compoundSet);
    return new ProteinSequence(sequence, compoundSet);
  }

  public AbstractSequence<AminoAcidCompound> getSequence(
      SequenceProxyLoader<AminoAcidCompound> proxyLoader, long index) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  public AbstractSequence<AminoAcidCompound> getSequence(
      List<AminoAcidCompound> list) {
    throw new UnsupportedOperationException("Not supported yet.");
  }
}
