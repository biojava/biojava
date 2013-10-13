/**
 * 
 */
package org.biojava3.core.sequence.io;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Collection;

import org.biojava3.core.sequence.io.template.GenbankHeaderFormatInterface;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.util.StringManipulationHelper;


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
