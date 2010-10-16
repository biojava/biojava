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
 * A no-frills implementation of the OrthoPair interface
 *
 * @author David Huen
 */
public class SimpleOrthoPair implements OrthoPair
{
    private Orthologue first;
    private Orthologue second;
    private SimilarityType type;
    private double percentIdentity = 0.0;
    private String ref = null;

    /**
     * constructor for the computed form
     * of an orthology relationship.
     *     
     */
    public SimpleOrthoPair(
        Orthologue first, 
        Orthologue second, 
        SimilarityType type, 
        double percentIdentity
        )
    {
        // validate the parameters
        if ((first == null) || (second == null)) throw new IllegalArgumentException();
        if ((type != SimilarityType.MULTIPLE) && (type != SimilarityType.TWIN))
            throw new IllegalArgumentException();

        // we always store the orthologies in ascending Taxon ID
        if (first.getTaxonID() < second.getTaxonID()) {
            this.first = first;
            this.second = second;
        }
        else {
            this.first = second;
            this.second = first;
        }

        this.type = type;
        this.percentIdentity = percentIdentity;
    }

    /**
     * constructor for the curated form
     * of an orthology relationship
     */
    public SimpleOrthoPair(
        Orthologue first,
        Orthologue second,
        String ref
        )
    {
        // validate the parameters
        if ((first == null) || (second == null)) throw new IllegalArgumentException();
        if (ref == null) throw new IllegalArgumentException();

        this.first = first;
        this.second = second;
        this.type = SimilarityType.CURATED;
    }

    public Orthologue getFirstOrthologue() { return first; }
    public Orthologue getSecondOrthologue() { return second; }
    public SimilarityType getSimilarity() { return type; }
    public double getPercentIdentity() { return percentIdentity; }
    public String getRef() { return ref; }

    public boolean equals(Object o)
    {
        if (!(o instanceof OrthoPair)) return false;

        OrthoPair other = (OrthoPair) o;

        // we do not need to reverse the relationship
        // since we store the entries in ascending taxonID
        // and Homologene does not model same-species paralogy.

        if (other.getFirstOrthologue().equals(first)) return false;
        if (other.getSecondOrthologue().equals(second)) return false;

        if (other.getSimilarity() != type) return false;

        // we do not check other fields as uniqueness ought to
        // have been determined already.  There should not
        // be another OrthoPair with the same orthologues
        // and SimilarityType.

        return true;
    }
}
