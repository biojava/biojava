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
 * Created on May 22, 2007
 * 
 */

package org.biojava.dasobert.feature;

import java.awt.Color;

public class SegmentImpl extends AbstractSegment  {

    public SegmentImpl() {
    	super();
        start = 0 ;
        end   = 0 ;
        name  = "Unknown";
        color = Color.white ;
        txtColor = "white" ;
        parent = null ;
        note = "";
    }
    
    public boolean equals(Segment s){
        if ( s == null)
            return false;
        
        if (    ( start == s.getStart() ) &&
                ( end == s.getEnd() ) &&
                ( name.equals(s.getName()))               
                )      
        {
            if ( note == null) {
                if (s.getNote() == null)
                    return true;
            } else {
                if (s.getNote() != null){
                    if (s.getNote().equals(note))
                        return true;
                }
            }
            
        }
        
        
        return false;
    }
    
    public Object clone(){
        
        Segment s = new SegmentImpl();
        s.setStart(start);
        s.setEnd(end);
        s.setName(name);
        s.setColor(color);
        s.setTxtColor(txtColor);
        s.setNote(note);
        return s;
        
    }
}
