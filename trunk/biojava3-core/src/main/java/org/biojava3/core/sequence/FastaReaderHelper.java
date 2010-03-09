/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.List;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaReaderHelper {

  /**
   * Read a fasta file containing amino acids with setup that would handle most
   * cases.
   *
   * @param file
   * @return
   * @throws Exception
   */

  public static List<ProteinSequence> readFastaProteinSequence(File file)
      throws Exception {
    FileInputStream inStream = new FileInputStream(file);
    List<ProteinSequence> proteinSequences = readFastaProteinSequence(inStream);
    inStream.close();
    return proteinSequences;
  }

  /**
   * Read a fasta file containing amino acids with setup that would handle most
   * cases. User is responsible for closing InputStream because you opened it
   *
   * @param inStream
   * @return
   * @throws Exception
   */

  public static List<ProteinSequence> readFastaProteinSequence(
      InputStream inStream) throws Exception {
    FastaReader<ProteinSequence,AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(
        inStream,
        new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(),
        new ProteinSequenceCreator(AminoAcidCompoundSet
            .getAminoAcidCompoundSet()));
    return fastaReader.process();
  }

}
