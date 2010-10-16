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
 *    WriterOutputHandler.java
 */
package org.biojava.utils.process;

import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Writer;


/**
 * Simple {@linkplain org.biojava.utils.process.OutputHandler output handler} 
 * that writes the output of an external process to an 
 * {@linkplain org.biojava.utils.process.ReaderWriterPipe#getWriter() writer}.
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class WriterOutputHandler extends ReaderWriterPipe 
implements OutputHandler {
    
    /* PRIVATE FIELDS */

    /**
     * The input stream for the external process.
     */
    private InputStream input;

    /**
     * Initializes the writer output handler.
     * @param writer the writer to which to write the output of the external
     * process. May be <code>null</code>.
     * @param tag a tag for logging. May be <code>null</code>.
     */
    public WriterOutputHandler(Writer writer, String tag) {
        super(null, writer, tag);
    }
    
    /* INTERFACE OutputHandler */

    /**
     * {@inheritDoc}
     */
    public void setInput(InputStream input) {
        this.input = input;
        if (input != null) {
            setReader(new InputStreamReader(input));
        }
    }

    /**
     * {@inheritDoc}
     */
    public InputStream getInput() {
        return input;
    }

}
