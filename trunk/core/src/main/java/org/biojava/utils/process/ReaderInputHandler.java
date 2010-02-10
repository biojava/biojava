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
 *    ReaderInputHandler.java
 */
package org.biojava.utils.process;

import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.util.logging.Logger;


/**
 * Reader {@linkplain org.biojava.utils.process.InputHandler input handler} 
 * that reads the input for an external process from a 
 * {@linkplain org.biojava.utils.process.ReaderWriterPipe#getReader() reader}. The
 * {@linkplain org.biojava.utils.process.StreamPipe#getOutput() output stream} for
 * the input of the external process is closed after the reader is read
 * to its end.
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class ReaderInputHandler extends ReaderWriterPipe implements
        InputHandler {
    
    /* STATIC FIELDS */

    /**
     * The class logger.
     */
    private static final Logger LOGGER = Logger
            .getLogger(ReaderInputHandler.class.getName());
    
    /* PRIVATE FIELDS */

    /**
     * The output stream for the external process.
     */
    private OutputStream output = null;
    
    /* PUBLIC CONSTRUCTORS */

    /**
     * Initializes the reader input handler.
     * @param reader the reader from which to read the input for the external
     * process. May be <code>null</code>.
     * @param tag a tag for logging. May be <code>null</code>. 
     */
    public ReaderInputHandler(Reader reader, String tag) {
        super(reader, null, tag);
    }
    
    /* INTERFACE OutputHandler */

    /**
     * {@inheritDoc}
     */
    public void setOutput(OutputStream output) {
        this.output = output;
        if (output != null) {
            setWriter(new OutputStreamWriter(output));
        }
    }
    
    /* INTERFACE OutputHandler */

    /**
     * {@inheritDoc}
     */
    public OutputStream getOutput() {
        return output;
    }
    
    /* INTERFACE Runnable */

    /**
     * {@inheritDoc}
     */
    public void run() {
        super.run();
        try {
            getWriter().close();
        } catch (IOException e) {
            LOGGER.severe(e.toString());
        }
    }
}
