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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import org.biojava3.core.io.ThingReader;
import org.biojava3.core.io.ThingReceiver;

/**
 * Can supply data.
 * @author Richard Holland
 * @since 3.0
 */
public abstract class FASTAReader implements ThingReader {

    private FASTAReceiver fastaReceiver;
    private String nextLine = null;
    private boolean eof = false;

    /**
     * Subclasses implement this to return the next line of data.
     * Should return null if EOF is hit.
     * @return the next line.
     * @throws IOException if the line could not be read.
     */
    protected abstract String readLine() throws IOException;

    public boolean canReadNextThing() throws IOException {
        if (this.nextLine == null && !this.eof) {
            do {
                this.nextLine = this.readLine();
                this.eof = this.nextLine == null;
            } while (!this.eof && !this.isHeaderLine(this.nextLine));
        }
        return !this.eof;
    }

    /**
     * Checks if the line is a header line.
     * @param line the line to check.
     * @return true if it is.
     */
    private boolean isHeaderLine(String line) {
        return line.startsWith(">");
    }

    public void readNextThing() throws IOException {
        // This if-check also moves our line pointer to the
        // next description line, if it's not already there.
        if (!this.canReadNextThing()) {
            throw new IOException("Attempt to parse past last item");
        }
        // Description line is current next line, from second char onwards.
        // (First char is the > symbol).
        this.fastaReceiver.setDescriptionLine(this.nextLine.substring(1));
        // Sequence is everything else.        
        this.fastaReceiver.setSequence(new StringBuilder()); // Initially empty.
        boolean allDone = false;
        do {
            this.nextLine = this.readLine();
            if (this.nextLine == null) {
                this.eof = true;
            } else {
                if (this.isHeaderLine(this.nextLine)) {
                    allDone = true;
                } else {
                    this.fastaReceiver.appendSequence(this.nextLine.trim());
                }
            }
        } while (!this.eof && !allDone);
    }

    public void setNextThingReceiver(ThingReceiver thingReceiver) {
        if (!(thingReceiver instanceof FASTAReceiver)) {
            throw new IllegalArgumentException("Cannot set a receiver which is not a FASTAReceiver");
        }
        this.fastaReceiver = (FASTAReceiver) thingReceiver;
    }
    
    /**
     * Reads FASTA from a Reader.
     */
    public static class FASTAReaderReader extends FASTAReader {
        private Reader r;
        private BufferedReader br;

        /** 
         * Read FASTA from the given Reader.
         * @param r the reader.
         */
        public FASTAReaderReader(Reader r) {
            this.r = r;
            this.br = new BufferedReader(this.r);
        }

        /**
         * Subclasses use this to get at the raw reader.
         * @return the reader.
         */
        protected Reader getReader() {
            return this.r;
        }
        
        @Override
        protected String readLine() throws IOException {
           return this.br.readLine();
        }
        
        public void close() throws IOException {
        }
    }
    
    /**
     * Reads FASTA from an InputStream.
     */
    public static class FASTAStreamReader extends FASTAReaderReader {
        public FASTAStreamReader(InputStream is) {
            super(new InputStreamReader(is));
        }
    }

    /**
     * Reads FASTA from a File.
     */
    public static class FASTAFileReader extends FASTAReaderReader {
        public FASTAFileReader(File f) throws IOException {
            super(new FileReader(f));
        }

        @Override
        public void close() throws IOException {
            super.close();
            this.getReader().close();
        }
    }
}
