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
 * Created on Jun 10, 2005
 *
 */
package org.biojava.dasobert.eventmodel;

/** An interface fore events related to selections of sequence
 * position, sequence range and locking of the selection.
 *
 * @author Andreas Prlic
 *
 */
public interface SequenceListener 
extends ObjectListener{
    
    /* select a certain sequence position */
    public void selectedSeqPosition(int position);
    
    /** select a certain range of a sequence 
     * @param start the start
     * @param end the end of the range
     * */
    public void selectedSeqRange(int start, int end);
    
    /** the current selecetion is locked and can not be changed 
     * @param flag true if selection should be locked
     * */
    public void selectionLocked(boolean flag);
    
    public void newSequence(SequenceEvent e);
    
    /** clear what has been selected
     * 
     *
     */
    public void clearSelection();
}
