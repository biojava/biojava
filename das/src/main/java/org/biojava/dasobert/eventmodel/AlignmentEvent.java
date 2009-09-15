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
 * Created on Nov 20, 2005
 *
 */
package org.biojava.dasobert.eventmodel;

import org.biojava.bio.program.das.dasalignment.Alignment;
import org.biojava.dasobert.dasregistry.Das1Source;

public class AlignmentEvent 
extends AbstractDasEvent{
	
	
    String accessionCode;
    Alignment alignment ;
    Alignment[] allAlignments;
    
    public AlignmentEvent(String accessionCode, 
    		Alignment ali,
    		Alignment[] allAlignments) {
        super();
        this.accessionCode = accessionCode;        
        alignment=ali;
        this.allAlignments = allAlignments;

    }
    
    /** create a new AlignmentEvent
     * 
     * @param accessionCode the accession code of the alignment
     * @param ali 
     * @param allAlignments 
     * @param source
     */
    public AlignmentEvent(String accessionCode, 
    		Alignment ali,
    		Alignment[] allAlignments,
    		Das1Source source) {
        super();
        this.accessionCode = accessionCode;        
        alignment=ali;
        this.allAlignments = allAlignments;
        this.dasSource = source;

    }
    
    public String getAccessionCode() {
        return accessionCode;
    }



    public void setAccessionCode(String accessionCode) {
        this.accessionCode = accessionCode;
    }



    /** get the first alignment 
     * @return Alignment the alignment nr 1*/
    public Alignment getAlignment(){
        return alignment;
    }
    
    /** get all alignments 
     * 
     * @return Alignment[] all alignments
     * */
    public Alignment[] getAllAlignments(){
        return allAlignments;
    }

}
