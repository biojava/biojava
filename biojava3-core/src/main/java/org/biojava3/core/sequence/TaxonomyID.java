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
 * Created on DATE
 *
 */
package org.biojava3.core.sequence;

/**
 * A sequence can be associated with a species or Taxonomy ID
 * @author Scooter Willis
 */
public class TaxonomyID {


    private String id = null;
    DataSource dataSource = DataSource.UNKNOWN;

    public TaxonomyID(String id, DataSource dataSource) {
        this.id = id;
        this.dataSource = dataSource;
    }

    /**
     * @return the id
     */
    public String getID() {
        return id;
    }

    /**
     * @return the source
     */
    public DataSource getDataSource() {
        return dataSource;
    }
}
