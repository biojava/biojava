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
package org.biojava3.file.fasta;

import java.io.Serializable;

import org.biojava3.core.io.ThingBuilder;
import org.biojava3.core.io.ThingEmitter;
import org.biojava3.core.io.ThingFormat;
import org.biojava3.core.io.ThingReader;
import org.biojava3.core.io.ThingWriter;
import org.biojava3.file.fasta.FASTAReader.FASTAFileReader;
import org.biojava3.file.fasta.FASTAReader.FASTAReaderReader;
import org.biojava3.file.fasta.FASTAReader.FASTAStreamReader;
import org.biojava3.file.fasta.FASTAWriter.FASTAFileWriter;
import org.biojava3.file.fasta.FASTAWriter.FASTAStreamWriter;
import org.biojava3.file.fasta.FASTAWriter.FASTAWriterWriter;

/**
 * Represents all you need to know about a FASTA file's contents.
 * 
 * @author Richard Holland
 * @since 3.0
 */
public interface FASTA extends Serializable {

	/**
	 * Use this as a parameter to ThingParserFactory to make it parse FASTA
	 * using the simpler method.
	 */
	public static final ThingFormat<FASTA> format = new ThingFormat<FASTA>() {

		public Class<? extends ThingBuilder<FASTA>> getBuilderClass() {
			return FASTABuilder.class;
		}

		public Class<? extends ThingEmitter<FASTA>> getEmitterClass() {
			return FASTAEmitter.class;
		}

		public Class<? extends ThingReader> getFileReaderClass() {
			return FASTAFileReader.class;
		}

		public Class<? extends ThingReader> getStreamReaderClass() {
			return FASTAStreamReader.class;
		}

		public Class<? extends ThingReader> getReaderReaderClass() {
			return FASTAReaderReader.class;
		}

		public Class<? extends ThingWriter> getFileWriterClass() {
			return FASTAFileWriter.class;
		}

		public Class<? extends ThingWriter> getStreamWriterClass() {
			return FASTAStreamWriter.class;
		}

		public Class<? extends ThingWriter> getWriterWriterClass() {
			return FASTAWriterWriter.class;
		}
	};

	/**
	 * Sets the description line on this FASTA object.
	 * 
	 * @param descLine
	 *            the description line.
	 */
	public void setDescriptionLine(String descLine);

	/**
	 * Gets the description line of this FASTA object.
	 * 
	 * @return the description line.
	 */
	public String getDescriptionLine();

	/**
	 * Sets the sequence for this FASTA object.
	 * 
	 * @param sequence
	 *            the sequence.
	 */
	public void setSequence(CharSequence seq);

	/**
	 * Gets the sequence for this FASTA object.
	 * 
	 * @return the sequence.
	 */
	public CharSequence getSequence();

	/**
	 * A simple implementation of the FASTA object type.
	 */
	public static class FASTAImpl implements FASTA {
		
		private static final long serialVersionUID = 1L;
		
		private String descLine = "";
		private CharSequence sequence = "";

		public String getDescriptionLine() {
			return this.descLine;
		}

		public CharSequence getSequence() {
			return this.sequence;
		}

		public void setDescriptionLine(String descLine) {
			if (descLine == null) {
				throw new NullPointerException(
						"Description line cannot be null.");
			}
			this.descLine = descLine;
		}

		public void setSequence(CharSequence seq) {
			if (seq == null) {
				throw new NullPointerException("Sequence cannot be null.");
			}
			this.sequence = seq;
		}
	}
}