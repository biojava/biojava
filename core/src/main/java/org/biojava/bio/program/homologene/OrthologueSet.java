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

package org.biojava.bio.program.homologene;

import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Interface for classes that store and manipulate
 * orthologues.
 * <p>
 * You cannot create Orthologues here, just
 * work with them.
 *
 * @author David Huen
 * * @author Matthew Pocock
 */
public interface OrthologueSet
{

    public static final ChangeType MODIFY = 
        new ChangeType("OrthologueSet modified",
            "org.biojava.bio.program.homologene.OrthologueSet",
            "MODIFY");

    /**
     * An iterator for the contents of
     * an OrthologueSet.
     */
    public interface Iterator
    {
        public boolean hasNext();
        public Orthologue nextOrthologue();
    }

    /*
     * Retrieve an orthologue from the set.
     */
    public Orthologue getOrthologue(String homologeneID);

    /**
     * Add an orthologue to the set.
     */
    public void addOrthologue(Orthologue ortho) throws ChangeVetoException;

    /**
     * Remove an orthologue from the set.
     */
    public void removeOrthologue(Orthologue ortho) throws ChangeVetoException;

    /**
     * Return an iterator to the contents of the set.
     */
    public Iterator iterator();

    /**
     * Filter the contents of a set.
     */
    public OrthologueSet filter(OrthologueFilter filter);
}
