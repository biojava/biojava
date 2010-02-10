/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of either the BSD licence or the GNU Lesser General
 * Public Licence.  These should be distributed with the code. 
 * If you do not have copies see:
 *
 *      http://www.opensource.org/licenses/bsd-license.php
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
 
package org.biojava.utils;

/**
 * Exception which is thrown if a child process managed by <code>ProcessTools</code>
 * exceeds a specified time limit.  Generally indicates that the results should
 * not be trusted!
 *
 * @author Thomas Down
 * @since 1.4
 */
 
public class ProcessTimeoutException extends Exception {
    private final int returnCode;
    
    public ProcessTimeoutException(int rc) {
        super();
        this.returnCode = rc;
    }
    
    public ProcessTimeoutException(int rc, String message) {
        super(message);
        this.returnCode = rc;
    }
    
    /**
     * Get the return code from the dying child process.
     */
    
    public int getReturnCode() {
        return returnCode;
    }
}