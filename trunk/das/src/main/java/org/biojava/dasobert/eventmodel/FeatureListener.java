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
 * Created on Oct 28, 2005
 *
 */
package org.biojava.dasobert.eventmodel;

//import org.biojava.spice.multipanel.eventmodel.FeatureEvent;

/**  a feature listener that returns the raw features as returned by a DAS source.
 * 
 */  
public interface FeatureListener {
    
    /** new features have been returned from the Annotation server 
     * 
     * @param e
     */
    public void newFeatures(FeatureEvent e);
    
    /** the server says that he is busy and we should try again in x seconds
     * 
     * @param e
     */
    public void comeBackLater(FeatureEvent e);
    
}

    
    
    

