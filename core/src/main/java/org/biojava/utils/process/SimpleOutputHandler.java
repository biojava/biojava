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
 *    SimpleOutputHandler.java
 */
package org.biojava.utils.process;

import java.io.OutputStream;


/**
 * Simple {@linkplain org.biojava.utils.process.OutputHandler output handler} 
 * that pipes the output of an external process to an 
 * {@linkplain org.biojava.utils.process.StreamPipe#getOutput() output stream}.
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class SimpleOutputHandler extends StreamPipe implements OutputHandler {

    /**
     * Initializes the simple output handler.
     * @param output the output stream to which to write the output of the 
     * external process. May be <code>null</code>.
     * @param tag a tag for logging. May be <code>null</code>.
     */
    public SimpleOutputHandler(OutputStream output, String tag) {
        super(null, output, tag);
    }

}
