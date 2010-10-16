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

/*
 *    ReaderWriterPipe.java
 */
package org.biojava.utils.process;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * A {@linkplain java.lang.Runnable multi threaded} class
 * which pipes the contents of an input reader to an output 
 * writer. 
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class ReaderWriterPipe implements Runnable {
    
    /* STATIC FIELDS */
    
    /**
     * The class logger.
     */
    private static final Logger LOGGER = 
        Logger.getLogger(ReaderWriterPipe.class.getName());
    
    /* PRIVATE FIELDS */

    /**
     * The reader from which to read.
     */
    private Reader reader;

    /**
     * The writer to which to write.
     */
    private Writer writer;

    /**
     * A tag for logging.
     */
    private String tag;
    
    /* PUBLIC CONSTRUCTORS */

    /**
     * Initializes the reader writer pipe.
     * @param reader the reader from which to read. May be <code>null</code>.
     * @param writer the writer to which to write. May be <code>null</code>.
     * @param tag a tag for loggging. May be <code>null</code>.
     */
    public ReaderWriterPipe(Reader reader, Writer writer, String tag) {
        setReader(reader);
        setWriter(writer);
        this.tag = tag;
    }
    
    /* PUBLIC PROPERTIES */

    /**
     * Gets the reader.
     * @return the reader from which to read. May be <code>null</code>.
     */
    public Reader getReader() {
        return reader;
    }

    /**
     * Gets the writer.
     * @return the writer to which to write. May be <code>null</code>.
     */
    public Writer getWriter() {
        return writer;
    }

    /**
     * Sets the reader.
     * @param reader the reader from which to read. May be <code>null</code>.
     */
    public void setReader(Reader reader) {
        this.reader = reader;
    }

    /**
     * Sets the writer.
     * @param writer the writer to which to write. May be <code>null</code>.
     */
    public void setWriter(Writer writer) {
        this.writer = writer;
    }

    /* INTERFACE Runnable */

    /**
     * {@inheritDoc}
     */
    public void run() {

        LOGGER.entering(getClass().getName(), "run");

        if (reader != null) {
            try {

                BufferedWriter bout = null;
                if (writer != null) {
                    bout = new BufferedWriter(writer);
                }
                BufferedReader bin = new BufferedReader(reader);
                boolean log = LOGGER.isLoggable(Level.FINEST);
                String line = null;
                while ((line = bin.readLine()) != null) {
                    if (bout != null) {
                        if (log) {
                            if (tag == null) {
                                LOGGER.finest(line);
                            } else {
                                LOGGER.finest("<" + tag + "> " + line);
                            }
                        }
                        bout.write(line);
                        bout.newLine();
                        bout.flush();
                    }
                }

            } catch (Exception e) {
                LOGGER.severe(e.toString());
            }
        }

        LOGGER.exiting(getClass().getName(), "run");
    }
}
