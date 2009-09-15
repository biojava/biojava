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
 *    OutputHandler.java
 */
package org.biojava.utils.process;

import java.io.InputStream;

/**
 * Interface to a {@linkplain java.lang.Runnable threadable} output handler for
 * an {@linkplain org.biojava.utils.process.ExternalProcess external process}. 
 * The output handler is used to collect the output of the 
 * {@linkplain java.lang.Process#getOutputStream()() STDOUT} output and/or the
 * {@linkplain java.lang.Process#getErrorStream() STDERR} of an external
 * process. 
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public interface OutputHandler extends Runnable {

    /**
     * Sets the input stream. The output of the external process is written to
     * this input stream.
     * @param input the input stream. May be <code>null</code>.
	 */
	void setInput(InputStream input);
    
    /**
     * Gets the input stream. The output of the external process is written to
     * this input stream.
     * @return the input stream. May be <code>null</code>.
     */
    /*@pure*/ InputStream getInput();
}
