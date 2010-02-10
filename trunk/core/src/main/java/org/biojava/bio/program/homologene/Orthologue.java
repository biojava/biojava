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
 * this entry contains data about the orthologue.
 */
public interface Orthologue
{
    /**
     * return the title used by Genbank for this protein
     */
    public String getTitle();

    /**
     * return the Taxon associated with this orthologue
     */
    public Taxon getTaxon();

    /**
     * a convenience method to return the TaxonID for thsi orhtologue.
     */
    public int getTaxonID();

    /**
     * get the locus ID associated with this orthologue.
     * It can be null.
     */
    public String getLocusID();

    /**
     * get the Homologene ID.  This is unique and
     * always defined.
     */
    public String getHomologeneID();

    /**
     * get the Accession ID associated with this orthologue.
     */
    public String getAccession();

    public void setTitle(String title);
}
