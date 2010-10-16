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
 * an interface for Homologene dataset Builders
 *
 * @author David Huen
 */
public interface HomologeneBuilder
{
    public String TAXONID = "TaxonID";
    public String LOCUSID = "LocusID";
    public String HOMOID = "HomologeneID";
    public String ACCESSION = "Accession";

    public String SIMILARITYTYPE = "SimilarityType";
    public String PERCENTIDENTITY = "PercentIdentity";
    public String REFERENCE = "Reference";

    public String TWIN = "twin";
    public String MULTIPLE = "multiple";
    public String CURATED = "curated";


    /**
     * indicates start of data for a HomologeneDB
     */
    public void startDB();

    /**
     * indicates start of data for a OrthoPairSet
     */
    public void startGroup();

    /**
     * indicates start of data for an OrthoPair
     */
    public void startOrthoPair();

    /**
     * indicates start of data for an orthologue
     */
    public void startOrthologue();

    /**
     * add a property to the current Orthologue
     */
    public void addOrthologueProperty(String key, String value);

    /**
     * end of data for this Orthologue
     */
    public void endOrthologue();
 
    /**
     * add a property to the current OrthoPair
     */
    public void addOrthoPairProperty(String key, String value);

    /**
     * end of data for this OrthoPair
     */
    public void endOrthoPair();

    /**
     * add title information to an Orthologue
     * (this is not in enclosed in the Orthologue element
     * because it comes completely separate in the Homologene
     * data files.  Go figger.)
     */
    public void addTitle(int taxonID, String homologeneID, String title);

    /**
     * end of data for group
     */
    public void endGroup();

    /**
     * end of data for DB
     */
    public void endDB();

    /**
     * retrieve the DB that has just been built
     */
    public HomologeneDB getDB();
}

