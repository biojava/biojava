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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;
import org.biojava3.core.io.ThingWriter;

/**
 * Receives data about FASTA and writes them to some form of output.
 * @author Richard Holland
 * @since 3.0
 */
public abstract class FASTAWriter implements ThingWriter, FASTAReceiver {

	private static final long serialVersionUID = 1L;
	
	private String descLine;
    private CharSequence seq;
    private PrintWriter writer;

    /**
     * Subclasses should use this to set up a print writer and call
     * {@link #setPrintWriter(java.io.PrintWriter)}.
     */
    public abstract void startGroupOfThings() throws IOException;

    /**
     * Subclasses sohuld use this to close the print writer.
     */
    public abstract void endGroupOfThings() throws IOException;

    /**
     * Subclasses should use this to pass in a writer for the format to
     * write to.
     * @param writer the writer.
     */
    protected void setPrintWriter(PrintWriter writer) {
        this.writer = writer;
    }

    /**
     * Get the print writer being used to write this file.
     * @return the print writer.
     */
    protected PrintWriter getPrintWriter() {
        return this.writer;
    }

    public String getDescriptionLine() {
        throw new UnsupportedOperationException();
    }

    public CharSequence getSequence() {
        throw new UnsupportedOperationException();
    }

    public void setDescriptionLine(String descLine) {
        this.descLine = descLine;
    }

    public void setSequence(CharSequence seq) {
        this.seq = seq;
    }
    
    public void appendSequence(CharSequence seq) {
        if (!(this.seq instanceof StringBuilder)) {
            this.seq = new StringBuilder(this.seq);
        }
        ((StringBuilder)this.seq).append(seq);
    }

    public void finishThing() throws IOException {
        this.writer.println(">" + this.descLine);
        for (int i = 0; i < this.seq.length(); i += 60) {
            this.writer.println(this.seq.subSequence(i, i + 60));
        }
    }

    public void startThing() throws IOException {
        this.descLine = "";
        this.seq = "";
    }

    /**
     * Write FASTA to a Writer. The writer will be flushed but
     * not closed on exit.
     */
    public static class FASTAWriterWriter extends FASTAWriter {

		private static final long serialVersionUID = 1L;
		
		private Writer w;

        /**
         * Set up a FASTA writer to a Writer.
         * @param w the Writer to write to.
         */
        public FASTAWriterWriter(Writer w) {
            this.w = w;
            this.setPrintWriter(new PrintWriter(this.w));
        }

        /**
         * Subclasses can get at the writer used in the constructor.
         * @return the writer.
         */
        protected Writer getWriter() {
            return this.w;
        }

        @Override
        public void startGroupOfThings() throws IOException {
        }

        @Override
        public void endGroupOfThings() throws IOException {
        }

        public void close() throws IOException {
            this.getPrintWriter().flush();
            this.getPrintWriter().close();
            this.w.flush();
        }
    }

    /**
     * Write FASTA to an OutputStream. The stream will be flushed but
     * not closed on exit.
     */
    public static class FASTAStreamWriter extends FASTAWriterWriter {

		private static final long serialVersionUID = 1L;
		
        /**
         * Set up a FASTA writer to an OutputStream.
         * @param os the OutputStream to write to.
         */
        public FASTAStreamWriter(OutputStream os) {
            super(new OutputStreamWriter(os));
        }
    }

    /**
     * Write FASTA to a File.
     */
    public static class FASTAFileWriter extends FASTAWriterWriter {

		private static final long serialVersionUID = 1L;
		
        /**
         * Set up a FASTA writer to a File.
         * @param file the File to write to.
         * @throws IOException if there were problems with the file.
         */
        public FASTAFileWriter(File file) throws IOException {
            super(new FileWriter(file));
        }

        @Override
        public void close() throws IOException {
            super.close();
            this.getWriter().close();
        }
    }
}