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

import org.biojava.bio.BioException;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/** 
 * <p>ARAlignment is an interface that defines methods for adding and
 * removing seqeunces from an Alignment.</p>
 *
 * @author David Waring
 */

public interface ARAlignment extends Alignment{
    
    public void addSequence(AlignmentElement ae) throws ChangeVetoException,BioException;
    public void removeSequence(Object label) throws ChangeVetoException;
    
    public static final ChangeType ADD_LABEL = new ChangeType(
        "Adding a sequence to the alignment",
        "org.biojava.bio.alignment.ARAlignment",
        "ADD_LABEL"
    );
    
    public static final ChangeType REMOVE_LABEL = new ChangeType(
        "Adding a sequence to the alignment",
        "org.biojava.bio.alignment.ARAlignment",
        "REMOVE_LABEL"
    );
    

}
