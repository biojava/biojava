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
 * Allows cross-references to other databases to be ranked.
 * @author Richard Holland
 * @author gwaldon
 * @see RankedCrossRefable
 * @see CrossRef
 * @since 1.5
 */
public interface RankedCrossRef extends Comparable,Changeable {
    
    public static final ChangeType RANK = new ChangeType(
            "This ranked crossreference's rank has changed",
            "org.biojavax.RankedCrossRef",
            "RANK"
            );
    
    /**
     * Return the cross reference associated with this object.
     * @return a crossref object.
     */
    public CrossRef getCrossRef();
    
    /**
     * Return the rank associated with the cross reference.
     * @return the rank.
     */
    public int getRank();
    
    /**
     * Set the rank associated with the cross reference.
     * @param rank the rank to use.
     * @throws ChangeVetoException if the new rank is unacceptable.
     */
    public void setRank(int rank) throws ChangeVetoException;
    
}
