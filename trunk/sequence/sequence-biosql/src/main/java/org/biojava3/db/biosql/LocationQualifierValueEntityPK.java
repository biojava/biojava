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
package org.biojava3.db.biosql;

import java.io.Serializable;
import javax.persistence.Column;
import javax.persistence.Embeddable;

/**
 * Embeddable PK for LocationQualifierValueEntity
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class LocationQualifierValueEntityPK implements Serializable {
    private static final long serialVersionUID = -7783560598176300100L;
    @Column(name = "LOCATION_ID", nullable = false)
    private int locationId;
    @Column(name = "TERM_ID", nullable = false)
    private int termId;

    public LocationQualifierValueEntityPK() {
    }

    public LocationQualifierValueEntityPK(int locationId, int termId) {
        this.locationId = locationId;
        this.termId = termId;
    }

    public int getLocationId() {
        return locationId;
    }

    public void setLocationId(int locationId) {
        this.locationId = locationId;
    }

    public int getTermId() {
        return termId;
    }

    public void setTermId(int termId) {
        this.termId = termId;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (int) locationId;
        hash += (int) termId;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof LocationQualifierValueEntityPK)) {
            return false;
        }
        LocationQualifierValueEntityPK other = (LocationQualifierValueEntityPK) object;
        if (this.locationId != other.locationId) {
            return false;
        }
        if (this.termId != other.termId) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.LocationQualifierValueEntityPK[locationId=" + locationId + ", termId=" + termId + "]";
    }

}
