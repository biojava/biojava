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
 *    StreamPipe.java
 */
package org.biojava.utils.process;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * A {@linkplain java.lang.Runnable multi threaded} class 
 * which pipes the contents of an input stream to an output stream. 
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class StreamPipe implements Runnable {
    
    /* STATIC FIELDS */
	
    /**
     * The class logger.
     */
    private static final Logger LOGGER = 
        Logger.getLogger(StreamPipe.class.getName());
    
    /* PRIVATE FIELDS */

    /**
     * The input stream from which to read.
     */
    private InputStream input = null;
    
	/**
	 * The output stream to which to write.
	 */
	private OutputStream output = null;
    
    /**
     * A tag for logging.
     */
    private String tag = null;
	
    /* PUBLIC CONSTRUCTORS */
    
	/**
     * Initializes the stream pipe.
	 * @param input the input stream from which to read. 
     * May be <code>null</code>.
	 * @param output the output stream to which to write
     * May be <code>null</code>.
	 * @param tag a tag which is used for logging the in- and output
     * May be <code>null</code>.
	 */
	public StreamPipe(InputStream input, OutputStream output, String tag) {
        setInput(input);
        setOutput(output);
		this.tag = tag;
	}
    
    /* PUBLIC PROPERTIES */

    /**
     * Gets the input stream
     * @return the input from which to read. May be <code>null</code>.
     */
    public InputStream getInput() {
        return input;
    }

    /**
     * Sets the input stream
     * @param input the input stream from which to read. May be 
     * <code>null</code>.
     */
    public void setInput(InputStream input) {
        this.input = input;
    }
    
    /**
     * Sets the output stream
     * @param output the output stream to which to write. May be 
     * <code>null</code>.
     */
    public void setOutput(OutputStream output) {
        this.output = output;
    }

    /**
     * Gets the output stream.
     * @return the output stream to which to write. May be <code>null</code>.
     */
    public OutputStream getOutput() {
        return output;
    }
    
    /* INTERFACE Runnable */

	/**
	 * {@inheritDoc}
	 */
	public void run() {
		
        LOGGER.entering(getClass().getName(), "run");
        
        if (input != null) {
      		try {  			
      			BufferedOutputStream bout = null;
      			if (output != null) {
                    bout = new BufferedOutputStream(output);
                }  
      			BufferedInputStream bin = new BufferedInputStream(input);
                boolean log = LOGGER.isLoggable(Level.FINEST);
      			byte[] buffer = new byte[1024];
      			int len;
      			while ((len = bin.read(buffer)) != -1) {
                    if (log) {
                        String data = new String(buffer, 0, len);
                        if (tag == null) {
                            LOGGER.finest(data);
                        } else {
                            LOGGER.finest("<" + tag + "> " + data);
                        }
                    }
                    if (bout != null) {
        				bout.write(buffer, 0, len);
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
