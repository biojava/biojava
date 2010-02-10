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


package org.biojava.bio.alignment;

import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <p>EditableAlignment is an interface that defines methods for
 * shifting bases within an Alignment.</p>
 *
 * <p>Shift should work as follows. Bases within a sequence can be
 * shifted to the right with offset > 1 to the left with offset <
 * 1. Shifting bases will be allowed if:</p>
 *
 * <li>1: Shift would remove only gaps on one side, they will be
 * replace with gaps on the other.
 *
 * <li>2: Shift is at the end of a sequence. It will add gaps if the
 * range location is less than the whole sequence.
 *
 * <li>Shifts that would delete bases will throw a
 * IllegalEditException
 *
 * <p>If the Alignment is an UnequalLengthAlignment it should be
 * acceptable to shift bases in such a way as to increase (or
 * decrease) the size of the overall length of the alignment,
 * i.e. shift them over the edge.</p>
 *
 * @author David Waring
 */
public interface EditableAlignment{
    
    /**
    * <P> edit() allows edits on an individual sequence, they should be reflected back
    * to the underlying SymbolList.
    */
    
    public void edit (Object label,Edit edit)throws ChangeVetoException;
    

    /**
    * loc in this case is the Alignment Location
    */
    public void shiftAtAlignmentLoc(Object label, Location loc, int offset) throws ChangeVetoException,IllegalAlignmentEditException,IndexOutOfBoundsException;
    
    
    /**
    * loc in this case is the SymbolList Location
    */
    public void shiftAtSequenceLoc(Object label, Location loc, int offset) throws ChangeVetoException,IllegalAlignmentEditException,IndexOutOfBoundsException;

    public static final ChangeType LOCATION = new ChangeType(
        "The location of a sequence is being changed",
        "org.biojava.bio.alignment.EditableAlignment",
        "LOCATION"
    );

    public static final ChangeType GAPS = new ChangeType(
        "The gap within a sequence are changing",
        "org.biojava.bio.alignment.EditableAlignment",
        "GAPS"
    );

 
}
