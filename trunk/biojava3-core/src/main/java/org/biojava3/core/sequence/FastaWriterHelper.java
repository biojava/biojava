/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.Collection;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.FastaWriter;
import org.biojava3.core.sequence.io.GenericFastaHeaderFormat;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaWriterHelper {

  /**
   * Write collection of protein sequences to a file
   *
   * @param file
   * @param proteinSequences
   * @throws Exception
   */

  public static void writeProteinSequence(File file,
      Collection<ProteinSequence> proteinSequences) throws Exception {
    FileOutputStream outputStream = new FileOutputStream(file);
    writeProteinSequence(outputStream, proteinSequences);
    outputStream.close();
  }

  public static void writeProteinSequence(OutputStream outputStream,
      Collection<ProteinSequence> proteinSequences) throws Exception {
    FastaWriter<ProteinSequence,AminoAcidCompound> fastaWriter = new FastaWriter<ProteinSequence,AminoAcidCompound>(
        outputStream, proteinSequences,
        new GenericFastaHeaderFormat<ProteinSequence,AminoAcidCompound>());
    fastaWriter.process();

  }
}
