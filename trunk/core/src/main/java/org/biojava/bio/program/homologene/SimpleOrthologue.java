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
public class SimpleOrthologue implements Orthologue
{
    private String title;
    private Taxon taxon;
    private String locusID;
    private String homologeneID;
    private String accession;

    public SimpleOrthologue(
        Taxon taxon,
        String locusID,
        String homologeneID,
        String accession
        )
    {
        setTaxon(taxon);
        setLocusID(locusID);
        setHomologeneID(homologeneID);
        setAccession(accession);
    }

    /**
     * this constructor does the Taxon lookup for you too
     */
    public SimpleOrthologue(
        int taxonID,
        String locusID,
        String homologeneID,
        String accession
        )
        throws IllegalArgumentException
    {
        // get corresponding Taxon
        taxon = HomologeneTools.getTaxon(taxonID);
        if (taxon == null) throw new IllegalArgumentException("Taxon with ID of " + taxonID + " does not exist.");

        setTaxon(taxon);
        setLocusID(locusID);
        setHomologeneID(homologeneID);
        setAccession(accession);
    }

    public String getTitle() { return title; }
    public Taxon getTaxon() { return taxon; }
    public int getTaxonID() { return taxon.getTaxonID(); }
    public String getLocusID() { return locusID; }
    public String getHomologeneID() { return homologeneID; }
    public String getAccession() { return accession; }

    public void setTitle(String title) { this.title = title; }
    void setTaxon(Taxon taxon) { this.taxon = taxon; }
    void setLocusID(String locusID) { this.locusID = locusID.trim(); }
    void setHomologeneID(String homologeneID) { this.homologeneID = homologeneID.trim(); }
    void setAccession(String accession) { this.accession = accession.trim(); }

    public boolean equals(Object o)
    {
        if (!(o instanceof Orthologue)) return false;

        // two Orthologues are only equal if they have identical data
        Orthologue other = (Orthologue) o;

        if (other.getTaxon() != taxon) return false;
        if (other.getLocusID() != locusID) return false;
        if (other.getHomologeneID() != homologeneID) return false;
        if (other.getAccession() != accession) return false;

        return true;
    }
}
