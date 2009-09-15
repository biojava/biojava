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
import java.util.Collection;
import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.OneToMany;
import javax.persistence.Table;

/**
 * The persitent entity bean that backs a Biodatabase implementation.
 * @author Mark Schreiber
 * @since 1.0
 */
@Entity
@Table(name = "BIODATABASE")
@NamedQueries({
    @NamedQuery(name = "BiodatabaseEntity.findByBiodatabaseEntityId", query = "SELECT b FROM BiodatabaseEntity b WHERE b.biodatabaseId = :biodatabaseId"),
    @NamedQuery(name = "BiodatabaseEntity.findByName", query = "SELECT b FROM BiodatabaseEntity b WHERE b.name = :name"), 
    @NamedQuery(name = "BiodatabaseEntity.findByAuthority", query = "SELECT b FROM BiodatabaseEntity b WHERE b.authority = :authority"), 
    @NamedQuery(name = "BiodatabaseEntity.findByDescription", query = "SELECT b FROM BiodatabaseEntity b WHERE b.description = :description")})
public class BiodatabaseEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "BIODATABASE_ID", nullable = false)
    private Integer biodatabaseId;
    @Column(name = "NAME", nullable = false)
    private String name;
    @Column(name = "AUTHORITY")
    private String authority;
    @Column(name = "DESCRIPTION")
    private String description;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "biodatabaseId")
    private Collection<BioentryEntity> bioentryCollection;

    public BiodatabaseEntity() {
    }

    public BiodatabaseEntity(Integer biodatabaseId) {
        this.biodatabaseId = biodatabaseId;
    }

    public BiodatabaseEntity(Integer biodatabaseId, String name) {
        this.biodatabaseId = biodatabaseId;
        this.name = name;
    }

    public Integer getBiodatabaseEntityId() {
        return biodatabaseId;
    }

    public void setBiodatabaseEntityId(Integer biodatabaseId) {
        this.biodatabaseId = biodatabaseId;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getAuthority() {
        return authority;
    }

    public void setAuthority(String authority) {
        this.authority = authority;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public Collection<BioentryEntity> getBioentryCollection() {
        return bioentryCollection;
    }

    public void setBioentryCollection(Collection<BioentryEntity> bioentryCollection) {
        this.bioentryCollection = bioentryCollection;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (biodatabaseId != null ? biodatabaseId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BiodatabaseEntity)) {
            return false;
        }
        BiodatabaseEntity other = (BiodatabaseEntity) object;
        if ((this.biodatabaseId == null && other.biodatabaseId != null) || (this.biodatabaseId != null && !this.biodatabaseId.equals(other.biodatabaseId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BiodatabaseEntity[biodatabaseId=" + biodatabaseId + "]";
    }

}
