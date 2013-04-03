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
 * Created on Aug 23, 2006
 *
 */
package org.biojava.bio.structure.gui.util;


/** a class that manages the conversion of sequence coordinate system to 
 * JPanel drawing coordinates
 * 
 * @author Andreas Prlic
 * @since 1.7
 * @version %I% %G%
 */
public class CoordManager {

    
    float scale;
    int chainLength;
    
    public CoordManager() {
        super();
        scale  = 1.0f;
        chainLength = 0;
    }
    
    public void setLength(int length){
        chainLength = length;
    }
    
    public void setScale(float scale){
        this.scale = scale;
    }

    
    /** start counting at 0...
     * 
     * @param panelPos
     * @return the sequence position
     */
    protected int getSeqPos(int panelPos){
     
        
        int seqPos = Math.round((panelPos - SequenceScalePanel.DEFAULT_X_START) / scale) ;
        if ( seqPos < 0)
            seqPos = 0;
        //int length = chainLength;
        //if ( seqPos >= length)
         //   seqPos = length-1;
        return seqPos;
    }
    
    protected int getPanelPos(int seqPos){
       
        if ( seqPos < 0 )
            seqPos = 0;
        
        //if ( seqPos >= length)
         //   seqPos = length-1;

        int aminosize = Math.round(1*scale);
        if ( aminosize < 1)
            aminosize = 1;
        
        int panelPos = Math.round(seqPos * scale) + SequenceScalePanel.DEFAULT_X_START ;
        return panelPos;
    }
    
    
}
