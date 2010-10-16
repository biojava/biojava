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

import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;


/**
 * A simple ranked comment designed to be used for BioEntry comments
 * in BioSQL. The Comment field is intended to be set by the constructor
 * and is immuntable; the Rank field is changeable.
 * @author Richard Holland
 * @author gwaldon
 * @see org.biojavax.bio.BioEntry
 * @since 1.5
 */
public interface Comment extends Comparable, Changeable {
    
    public static final ChangeType RANK = new ChangeType(
            "This comment's rank has changed",
            "org.biojavax.Comment",
            "RANK"
            );
    
    /**
     * Returns the comment part of this comment.
     * @return a comment.
     */
    public String getComment();
    
    /**
     * Returns the rank of this comment.
     * @return the rank.
     */
    public int getRank();
    
    /**
     * Sets the rank of this comment.
     * @param rank the rank to use.
     * @throws ChangeVetoException if the new rank is unacceptable.
     */
    public void setRank(int rank) throws ChangeVetoException;
    
}
