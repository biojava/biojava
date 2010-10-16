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


/**
 * Each HomologeneEntry represents a single
 * Homologene record that relates two
 * presumptive orthologues.
 */
public interface OrthoPair
{
    /**
     * gets the first orthologue in the orthology
     * relationship.
     * <p>
     * I agree that this access route is somewhat
     * sucky for a symmetric relationship 
     * but I'd hate returning a Set of
     * orthologues more and all the fiddling
     * would require.
     */
    public Orthologue getFirstOrthologue();

    /**
     * gets the first orthologue in the orthology
     * relationship.  This will be the one with the
     * lower TaxonID.
     */
    public Orthologue getSecondOrthologue();

    /**
     * gets the second orthologue in the
     * orthology relationship.  This will be
     * the one with the higher TaxonID.
     */
    public SimilarityType getSimilarity();

    /**
     * get percentage identity.
     */
    public double getPercentIdentity();

    /**
     * get reference to evidence for 
     * orthology.
     */
    public String getRef();
}
