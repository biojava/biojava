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

package org.biojavax;

/**
 * An exception that indicates that an attempt to resolve a <code>CrossRef</code>
 * has failed.
 *
 * @author Mark Schreiber
 * @since 1.5
 */
public class CrossReferenceResolutionException extends Exception{
    
    /** Creates a new instance of CrossReferenceResolutionException */
    public CrossReferenceResolutionException() {
        super();
    }
    
    /** 
     * Creates a new instance of CrossReferenceResolutionException with a
     * message.
     * @param message a description or reason for the exception.
     */
    public CrossReferenceResolutionException(String message) {
        super(message);
    }
    
    /** 
     * Creates a new instance of CrossReferenceResolutionException with a
     * message and a cause.
     * @param message a description or reason for the exception.
     * @param cause the exception that caused this one.
     */
    public CrossReferenceResolutionException(String message, Throwable cause){
        super(message, cause);
    }
    
    /** 
     * Creates a new instance of CrossReferenceResolutionException with a
     * cause.
     * @param cause the exception that caused this one.
     */
    public CrossReferenceResolutionException(Throwable cause) {
        super(cause);
    }
}
