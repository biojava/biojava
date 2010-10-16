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
 *    SimpleInputHandler.java
 */
package org.biojava.utils.process;

import java.io.IOException;
import java.io.InputStream;
import java.util.logging.Logger;


/**
 * Simple {@linkplain org.biojava.utils.process.InputHandler input handler} 
 * that reads the input for an external process from an 
 * {@linkplain org.biojava.utils.process.StreamPipe#getInput() input stream}. The
 * {@linkplain org.biojava.utils.process.StreamPipe#getOutput() output stream} for
 * the input of the external process is closed after the input stream is read
 * to its end.
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public final class SimpleInputHandler extends StreamPipe 
implements InputHandler {
    
    /* STATIC FIELDS */
    
    /**
     * The class logger.
     */
    private static final Logger LOGGER = 
        Logger.getLogger(SimpleInputHandler.class.getName());

    /* PUBLIC CONSTRUCTORS */
    
    /**
     * Initializes the simple input handler.
     * @param input the input stream from which to read the input for the 
     * external process. May be <code>null</code>.
     * @param tag a tag for logging. May be <code>null</code>.
     */
    public SimpleInputHandler(InputStream input, String tag) {
        super(input, null, tag);
    }
    
    /* INTERFACE Runnable */

    /**
     * {@inheritDoc}
     */
    public void run() {
        super.run();
        try {
            getOutput().close();
        } catch (IOException e) {
            LOGGER.severe(e.toString());
        }
    }

}
