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
import javax.persistence.EmbeddedId;
import javax.persistence.Entity;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.Table;

/**
 * Entity for LocationQualifierValueEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "LOCATION_QUALIFIER_VALUE")
@NamedQueries({@NamedQuery(name = "LocationQualifierValueEntity.findByLocationId", query = "SELECT l FROM LocationQualifierValueEntity l WHERE l.locationQualifierValuePK.locationId = :locationId"), @NamedQuery(name = "LocationQualifierValueEntity.findByTermId", query = "SELECT l FROM LocationQualifierValueEntity l WHERE l.locationQualifierValuePK.termId = :termId"), @NamedQuery(name = "LocationQualifierValueEntity.findByValue", query = "SELECT l FROM LocationQualifierValueEntity l WHERE l.value = :value"), @NamedQuery(name = "LocationQualifierValueEntity.findByIntValue", query = "SELECT l FROM LocationQualifierValueEntity l WHERE l.intValue = :intValue")})
public class LocationQualifierValueEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    protected LocationQualifierValueEntityPK locationQualifierValuePK;
    @Column(name = "VALUE", nullable = false)
    private String value;
    @Column(name = "INT_VALUE")
    private Integer intValue;
    @JoinColumn(name = "LOCATION_ID", referencedColumnName = "LOCATION_ID", insertable = false, updatable = false)
    @ManyToOne
    private LocationEntity location;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID", insertable = false, updatable = false)
    @ManyToOne
    private TermEntity term;

    public LocationQualifierValueEntity() {
    }

    public LocationQualifierValueEntity(LocationQualifierValueEntityPK locationQualifierValuePK) {
        this.locationQualifierValuePK = locationQualifierValuePK;
    }

    public LocationQualifierValueEntity(LocationQualifierValueEntityPK locationQualifierValuePK, String value) {
        this.locationQualifierValuePK = locationQualifierValuePK;
        this.value = value;
    }

    public LocationQualifierValueEntity(int locationId, int termId) {
        this.locationQualifierValuePK = new LocationQualifierValueEntityPK(locationId, termId);
    }

    public LocationQualifierValueEntityPK getLocationQualifierValuePK() {
        return locationQualifierValuePK;
    }

    public void setLocationQualifierValuePK(LocationQualifierValueEntityPK locationQualifierValuePK) {
        this.locationQualifierValuePK = locationQualifierValuePK;
    }

    public String getValue() {
        return value;
    }

    public void setValue(String value) {
        this.value = value;
    }

    public Integer getIntValue() {
        return intValue;
    }

    public void setIntValue(Integer intValue) {
        this.intValue = intValue;
    }

    public LocationEntity getLocation() {
        return location;
    }

    public void setLocation(LocationEntity location) {
        this.location = location;
    }

    public TermEntity getTerm() {
        return term;
    }

    public void setTerm(TermEntity term) {
        this.term = term;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (locationQualifierValuePK != null ? locationQualifierValuePK.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof LocationQualifierValueEntity)) {
            return false;
        }
        LocationQualifierValueEntity other = (LocationQualifierValueEntity) object;
        if ((this.locationQualifierValuePK == null && other.locationQualifierValuePK != null) || (this.locationQualifierValuePK != null && !this.locationQualifierValuePK.equals(other.locationQualifierValuePK))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.LocationQualifierValueEntity[locationQualifierValuePK=" + locationQualifierValuePK + "]";
    }

}
