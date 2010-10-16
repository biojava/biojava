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
public interface SimilarityType
{
    public final static SimilarityType MULTIPLE
        = new PlaceHolder("reciprocal best match between 3 or more organisms");
    public final static SimilarityType TWIN
        = new PlaceHolder("reciprocal best match between 2 organisms");
    public final static SimilarityType CURATED
        = new PlaceHolder("curated homology relationship");

    public class PlaceHolder implements SimilarityType
    {
        String description;

        private PlaceHolder(String description)
        {
            this.description = description;
        }

        public String getDescription()
        {
            return description;
        }
    }
}
