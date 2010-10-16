/*
 *                  BioJava development code
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
 * Created on Nov 1, 2005
 *
 */
package org.biojava.dasobert.eventmodel;

/** an interface for the listeners of new PDB code requested / new Uniprot code requested
 * 
 * @author Andreas Prlic
 *
 */
public interface ObjectListener {
    
    /** a new object has been requested 
     * 
     * @param accessionCode
     */
    public void newObjectRequested(String accessionCode);
    
    /** no object with that accessionCode has been found
     * 
     * @param accessionCode
     */
    public void noObjectFound(String accessionCode);
    
    
   // public void exceptionOccured(Exception e);
   
}
