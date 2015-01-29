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

import org.biojava3.core.util.Equals;
import org.biojava3.core.util.Hashcoder;

/**
 * Used in Sequences as the unique indentifier. If possible, set the {@link DataSource} to know the
 * source of the id. This allows a SequenceProxy to gather features or related sequences
 * Protein->Gene as an example. When parsing a Blast file it is also possible
 * to identify the type of ID
 *
 * @author Scooter Willis
 */
public class AccessionID {

    private String id = null;
    private DataSource source = DataSource.LOCAL;

    /**
     *
     */

    public AccessionID(){
        id = "";
        
    }

    /**
     *
     * @param id
     */
    public AccessionID(String id) {
        this.id = id.trim();
        this.source = DataSource.LOCAL;
    }

    /**
     *
     * @param id
     * @param source
     */
    public AccessionID(String id, DataSource source) {
        this.id = id.trim();
        this.source = source;
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
        return source;
    }

    @Override
    public boolean equals(Object o) {
        boolean equals = false;
        if (Equals.classEqual(this, o)) {
            AccessionID l = (AccessionID) o;
            equals = (Equals.equal(getID(), l.getID())
                    && Equals.equal(getDataSource(), l.getDataSource()));
        }
        return equals;
    }

    @Override
    public int hashCode() {
        int r = Hashcoder.SEED;
        r = Hashcoder.hash(r, getID());
        r = Hashcoder.hash(r, getDataSource());
        return r;
    }

 //   public void setDataSource(DataSource dataSource){
 //       source = dataSource;
 //   }

    @Override
    public String toString() {
        return id;
    }
}
