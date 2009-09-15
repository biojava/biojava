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

import org.biojava.utils.ChangeVetoException;

/**
 * Homologene is a NCBI dataset that curates sets
 * of orthologues from the reference model organisms.
 * <p>
 * This class is a Collection of methods for handling
 * data from the Homologene dataset.
 *
 * @author David Huen
 * @author Matthew Pocock
 */
public interface HomologeneDB
{
    /**
     * Create an orthologue.
     */
    public Orthologue createOrthologue(Taxon taxon, String locusID, String homologeneID, String accession)
        throws ChangeVetoException;

    /**
     * Create an orthologue.
     */
    public Orthologue createOrthologue(int taxonID, String locusID, String homologeneID, String accession)
        throws ChangeVetoException;

    /**
     * Returns an orthologue of specified ID.
     */
    public Orthologue getOrthologue(String homologeneID);

    /**
     * Create a computed orthology entry.
     */
    public OrthoPair createOrthoPair(Orthologue first, Orthologue second, SimilarityType type, double percentIdentity);

    /**
     * Create a curated orthology entry.
     */
    public OrthoPair createOrthoPair(Orthologue first, Orthologue second, String ref);

    /**
     * Create a Homologene Group.
     */
    public OrthoPairSet createOrthoPairSet();    

    /**
     * Get the HomologeneGroups in this database.
     */
    public OrthoPairCollection getOrthoPairSets();

    /**
     * Filter the database for a specified group.
     */
    public OrthoPairCollection filter(OrthoPairSetFilter filters);
}

