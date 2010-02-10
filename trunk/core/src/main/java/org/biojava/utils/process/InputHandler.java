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
 *    InputHandler.java
 */
package org.biojava.utils.process;

import java.io.OutputStream;

/**
 * Interface to a {@linkplain java.lang.Runnable threadable} input handler for
 * an {@linkplain org.biojava.utils.process.ExternalProcess external process}. 
 * The input handler is used to write the
 * {@linkplain java.lang.Process#getInputStream() STDIN} input of an
 * external process. 
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public interface InputHandler extends Runnable {

	/**
     * Sets the output stream. The input for the external process is read from
     * this output stream.
	 * @param output the output stream. May be <code>null</code>.
	 */
	void setOutput(OutputStream output);
	
    /**
     * Gets the output stream. The input for the external process is read from
     * this output stream.
     * @return the output stream. May be <code>null</code>.
     */
    /*@pure@*/ OutputStream getOutput();
}
