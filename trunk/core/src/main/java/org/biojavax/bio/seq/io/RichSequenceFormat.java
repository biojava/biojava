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

package org.biojavax.bio.seq.io;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojavax.Namespace;
import org.biojavax.bio.seq.RichSequence;


/**
 * Allows a file format to be read/written as RichSequences.
 * @author Richard Holland
 * @since 1.5
 */
public interface RichSequenceFormat extends SequenceFormat {
    /**
     * Check to see if a given file is in our format. Some formats may be
     * able to determine this by filename, whilst others may have to open the
     * file and read it to see what format it is in.
     * @param file  the <code>File</code> to check.
     * @return true if the file is readable by this format, false if not.
     * @throws IOException in case the file is inaccessible.
     */
    public boolean canRead(File file) throws IOException;
    
    /**
     * On the assumption that the file is readable by this format (not checked),
     * attempt to guess which symbol tokenization we should use to read it.
     * For formats that only accept one tokenization, just return it without
     * checking the file. For formats that accept multiple tokenizations, its
     * up to you how you do it.
     * @param file  the <code>File</code> object to guess the format of.
     * @return a <code>SymbolTokenization</code> to read the file with.
     * @throws IOException if the file is unrecognisable or inaccessible.
     */
    public SymbolTokenization guessSymbolTokenization(File file) throws IOException;
    
    /**
     * Check to see if a given stream is in our format. 
     * @param stream the <code>BufferedInputStream</code> to check.
     * @return true if the stream is readable by this format, false if not.
     * @throws IOException in case the stream is inaccessible.
     */
    public boolean canRead(BufferedInputStream stream) throws IOException;
    
    /**
     * On the assumption that the stream is readable by this format (not checked),
     * attempt to guess which symbol tokenization we should use to read it.
     * For formats that only accept one tokenization, just return it without
     * checking the stream. For formats that accept multiple tokenizations, its
     * up to you how you do it.
     * @param stream the <code>BufferedInputStream</code> object to guess the format of.
     * @return a <code>SymbolTokenization</code> to read the stream with.
     * @throws IOException if the stream is unrecognisable or inaccessible.
     */
    public SymbolTokenization guessSymbolTokenization(BufferedInputStream stream) throws IOException;
    
    /**
     * Sets the stream to write to.
     * @param os the PrintStream to write to.
     * @throws IOException if writing fails.
     */
    public void setPrintStream(PrintStream os);
    
    /**
     * Gets the print stream currently being written to.
     * @return the current print stream.
     */
    public PrintStream getPrintStream();
    
    /**
     * Informs the writer that we want to start writing. This will do any initialisation
     * required, such as writing the opening tags of an XML file that groups sequences together.
     * @throws IOException if writing fails.
     */
    public void beginWriting() throws IOException;
    
    /**
     * Informs the writer that are done writing. This will do any finalisation
     * required, such as writing the closing tags of an XML file that groups sequences together.
     * @throws IOException if writing fails.
     */
    public void finishWriting() throws IOException;
    
    /**
     * Reads a sequence from the given buffered reader using the given tokenizer to parse
     * sequence symbols. Events are passed to the listener, and the namespace used
     * for sequences read is the one given. If the namespace is null, then the default
     * namespace for the parser is used, which may depend on individual implementations
     * of this interface.
     * @param reader the input source
     * @param symParser the tokenizer which understands the sequence being read
     * @param listener the listener to send sequence events to
     * @param ns the namespace to read sequences into.
     * @return true if there is more to read after this, false otherwise.
     * @throws BioException in case of parsing errors.
     * @throws IllegalSymbolException if the tokenizer couldn't understand one of the
     * sequence symbols in the file.
     * @throws IOException if there was a read error.
     */
    public boolean readRichSequence(BufferedReader reader, SymbolTokenization symParser,
            RichSeqIOListener listener, Namespace ns) throws BioException, IllegalSymbolException, IOException;
    
    /**
     * Writes a sequence out to the outputstream given by beginWriting() using the default format of the
     * implementing class. If namespace is given, sequences will be written with that
     * namespace, otherwise they will be written with the default namespace of the
     * implementing class (which is usually the namespace of the sequence itself).
     * If you pass this method a sequence which is not a RichSequence, it will attempt to
     * convert it using RichSequence.Tools.enrich(). Obviously this is not going to guarantee
     * a perfect conversion, so it's better if you just use RichSequences to start with!
     * @param seq the sequence to write
     * @param ns the namespace to write it with
     * @throws IOException in case it couldn't write something
     */
    public void writeSequence(Sequence seq, Namespace ns) throws IOException;
    
    /**
     * Retrive the current line width. Defaults to 80.
     * @return the line width
     */
    public int getLineWidth();
    
    /**
     * Set the line width. When writing, the lines of sequence will never be longer than the line
     * width. Defaults to 80.
     * @param width the new line width
     */
    public void setLineWidth(int width);
    
    /**
     * Use this method to toggle reading of sequence data.
     * @param elideSymbols set to true if you <em>don't</em> want the sequence data.
     */
    public void setElideSymbols(boolean elideSymbols);
    
    /**
     * Is the format going to emit events when sequence data is read?
     * @return true if it is <em>not</em> otherwise false (false is default) .
     */
    public boolean getElideSymbols();
    
    /**
     * Use this method to toggle reading of feature data.
     * @param elideFeatures set to true if you <em>don't</em> want the feature data.
     */
    public void setElideFeatures(boolean elideFeatures);
    
    /**
     * Is the format going to emit events when feature data is read?
     * @return true if it is <em>not</em> otherwise false (false is default).
     */
    public boolean getElideFeatures();
    
    /**
     * Use this method to toggle reading of bibliographic reference data.
     * @param elideReferences set to true if you <em>don't</em> want the bibliographic reference data.
     */
    public void setElideReferences(boolean elideReferences);
    
    /**
     * Is the format going to emit events when bibliographic reference data is read?
     * @return true if it is <em>not</em> otherwise false (false is default) .
     */
    public boolean getElideReferences();
    
    /**
     * Use this method to toggle reading of comments data. Will also ignore remarks
     * lines in bibliographic references.
     * @param elideComments set to true if you <em>don't</em> want the comments data.
     */
    public void setElideComments(boolean elideComments);
    
    /**
     * Is the format going to emit events when comments data or remarks from
     * bibliographic references are read?
     * @return true if it is <em>not</em> otherwise false (false is default).
     */
    public boolean getElideComments();
    
    /**
     * Provides a basic format with simple things like line-widths precoded.
     */
    public abstract class BasicFormat implements RichSequenceFormat  {
        
        private int lineWidth = 80;
        private boolean elideSymbols = false;
        private boolean elideFeatures = false;
        private boolean elideComments = false;
        private boolean elideReferences = false;
        private PrintStream os;
        
        /**
         * {@inheritDoc}
         */
        public boolean canRead(File file) throws IOException {
            return false;
        }
        
        /**
         * {@inheritDoc}
         */
        public SymbolTokenization guessSymbolTokenization(File file) throws IOException {
            return RichSequence.IOTools.getDNAParser();
        }
        
        /**
         * {@inheritDoc}
         */
        public int getLineWidth() { return this.lineWidth; }
        
        /**
         * {@inheritDoc}
         */
        public void setLineWidth(int width) {
            if (width<1) throw new IllegalArgumentException("Width cannot be less than 1");
            this.lineWidth = width;
        }
        
        /**
         * {@inheritDoc}
         */
        public boolean getElideSymbols() { return this.elideSymbols; }
        
        /**
         * {@inheritDoc}
         */
        public void setElideSymbols(boolean elideSymbols) { this.elideSymbols = elideSymbols; }
        
        /**
         * {@inheritDoc}
         */
        public boolean getElideFeatures() { return this.elideFeatures; }
        
        /**
         * {@inheritDoc}
         */
        public void setElideFeatures(boolean elideFeatures) { this.elideFeatures = elideFeatures; }
        
        /**
         * {@inheritDoc}
         */
        public boolean getElideReferences() { return this.elideReferences; }
        
        /**
         * {@inheritDoc}
         */
        public void setElideReferences(boolean elideReferences) { this.elideReferences = elideReferences; }
        
        /**
         * {@inheritDoc}
         */
        public boolean getElideComments() { return this.elideComments; }
        
        /**
         * {@inheritDoc}
         */
        public void setElideComments(boolean elideComments) { this.elideComments = elideComments; }
        
        /**
         * {@inheritDoc}
         */
        public void setPrintStream(PrintStream os) {
            if (os==null) throw new IllegalArgumentException("Print stream cannot be null");
            this.os = os;
        }
        
        /**
         * {@inheritDoc}
         */
        public PrintStream getPrintStream() { return this.os; }
    }
    
    /**
     * Provides the basic implementation required for simple header/footer-less files such as Genbank.
     */
    public abstract class HeaderlessFormat extends BasicFormat {
        /**
         * {@inheritDoc}
         */
        public void beginWriting() throws IOException {}
        
        /**
         * {@inheritDoc}
         */
        public void finishWriting() throws IOException {}
    }
}
