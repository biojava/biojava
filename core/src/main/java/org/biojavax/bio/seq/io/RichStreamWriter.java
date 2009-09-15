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

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.SequenceIterator;
import org.biojavax.Namespace;


/**
 * Writes all of the sequences from a SequenceIterator to a stream with a
 * particular format.
 * This can be wired from a StreamReader to make a simple file-format conversion
 * utility, or can be used to write out the sequences in a database to disk.
 * @author Matthew Pocock
 * @author Richard Holland
 * @since 1.5
 */
public class RichStreamWriter {
    
    /**
     * The format to write in.
     */
    private RichSequenceFormat format;
    
    /**
     * The stream to write to.
     */
    private PrintStream os;
    
    /**
     * Write each of the sequences in ss to the stream in the given format.
     * @param ss  the SequenceIterator to loop over
     * @throws IOException if the stream has any problems
     */
    public void writeStream(SequenceIterator ss, Namespace ns)
    throws IOException {
        this.format.setPrintStream(this.os);
        this.format.beginWriting();
        while(ss.hasNext()) {
            try {
                this.format.writeSequence(ss.nextSequence(), ns);
            } catch (BioException se) {
                se.printStackTrace();
            }
        }
        this.format.finishWriting();
    }
    
    /**
     * Generate a new RichStreamWriter to the stream os and using format.
     * @param os  the OutputStream to write to
     * @param format the SequenceFormat to write with
     */
    public RichStreamWriter(OutputStream os, RichSequenceFormat format) {
        this.os = new PrintStream(os);
        this.format = format;
    }
}
