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
 */
/**
 *
 */
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.sequence.io.template.GenbankHeaderFormatInterface;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.util.StringManipulationHelper;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Collection;


/**
 * @author mckeee1
 *
 */
public class GenbankWriter<S extends Sequence<?>, C extends Compound> {
	int SEQUENCE_INDENT = 9;

	OutputStream os;
	Collection<S> sequences;
	GenbankHeaderFormatInterface<S, C> headerFormat;
	private int lineLength = 60;

	// byte[] lineSep = System.getProperty("line.separator").getBytes();
	/**
	 * Use default line length of 60
	 *
	 * @param os
	 * @param sequences
	 * @param headerFormat
	 */
	public GenbankWriter(OutputStream os, Collection<S> sequences,
			GenbankHeaderFormatInterface<S, C> headerFormat) {

		this.os = os;
		this.sequences = sequences;
		this.headerFormat = headerFormat;
	}

	/**
	 * Set custom lineLength
	 *
	 * @param os
	 * @param sequences
	 * @param headerFormat
	 * @param lineLength
	 */

	public GenbankWriter(OutputStream os, Collection<S> sequences,
			GenbankHeaderFormatInterface<S, C> headerFormat, int lineLength) {
		this.os = os;
		this.sequences = sequences;
		this.headerFormat = headerFormat;
		this.lineLength = lineLength;
	}

	/**
	 * Allow an override of operating system line separator for programs that
	 * needs a specific CRLF or CR or LF option
	 *
	 * @param lineSeparator
	 */

	public void process() throws Exception {
		// Loosely based on code from Howard Salis
		// TODO - Force lower case?
		// boolean closeit = false;
		PrintWriter writer = new PrintWriter(os);
		for (S sequence : sequences) {
			String header = headerFormat.getHeader(sequence);
			writer.format(header);
			writer.println();
			// os.write(lineSep);

			/*
			 * if isinstance(record.seq, UnknownSeq): #We have already recorded
			 * the length, and there is no need #to record a long sequence of
			 * NNNNNNN...NNN or whatever. if "contig" in record.annotations:
			 * self._write_contig(record) else: self.handle.write("ORIGIN\n")
			 * return
			 */

			String data = sequence.getSequenceAsString().toLowerCase();
			int seq_len = data.length();
			writer.println("ORIGIN");
			// os.write(lineSep);

			for (int line_number = 0; line_number < seq_len; line_number += lineLength) {
				writer.print(StringManipulationHelper.padLeft(
						Integer.toString(line_number + 1), SEQUENCE_INDENT));
				for (int words = line_number; words < Math.min(line_number
						+ lineLength, seq_len); words += 10) {
					if ((words + 10) > data.length()) {
						writer.print((" " + data.substring(words)));
					} else {
						writer.print((" " + data.substring(words, words + 10)));
					}
				}
				// os.write(lineSep);
				writer.println();
			}

			writer.println("//");

		}

		writer.flush();

	}

	/*
	 * public static void main(String[] args) { try { FileInputStream is = new
	 * FileInputStream("/Users/Scooter/scripps/dyadic/c1-454Scaffolds.faa");
	 *
	 *
	 * FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new
	 * FastaReader<ProteinSequence, AminoAcidCompound>(is, new
	 * GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(), new
	 * ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
	 * LinkedHashMap<String, ProteinSequence> proteinSequences =
	 * fastaReader.process(); is.close();
	 *
	 *
	 * // System.out.println(proteinSequences);
	 *
	 * FileOutputStream fileOutputStream = new
	 * FileOutputStream("/Users/Scooter/scripps/dyadic/c1-454Scaffolds_temp.faa"
	 * );
	 *
	 * BufferedOutputStream bo = new BufferedOutputStream(fileOutputStream);
	 * long start = System.currentTimeMillis(); FastaWriter<ProteinSequence,
	 * AminoAcidCompound> fastaWriter = new FastaWriter<ProteinSequence,
	 * AminoAcidCompound>(bo, proteinSequences.values(), new
	 * GenericFastaHeaderFormat<ProteinSequence, AminoAcidCompound>());
	 * fastaWriter.process(); bo.close(); long end = System.currentTimeMillis();
	 * System.out.println("Took " + (end - start) + " seconds");
	 *
	 * fileOutputStream.close();
	 *
	 *
	 * } catch (Exception e) { e.printStackTrace(); } }
	 */
	/**
	 * @return the lineLength
	 */
	public int getLineLength() {
		return lineLength;
	}

	/**
	 * @param lineLength
	 *            the lineLength to set
	 */
	public void setLineLength(int lineLength) {
		this.lineLength = lineLength;
	}

}
