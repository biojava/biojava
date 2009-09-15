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
 * Entity for DbxrefQualifierValue
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "DBXREF_QUALIFIER_VALUE")
@NamedQueries({@NamedQuery(name = "DbxrefQualifierValueEntity.findByDbxrefId", query = "SELECT d FROM DbxrefQualifierValueEntity d WHERE d.dbxrefQualifierValuePK.dbxrefId = :dbxrefId"), @NamedQuery(name = "DbxrefQualifierValueEntity.findByTermId", query = "SELECT d FROM DbxrefQualifierValueEntity d WHERE d.dbxrefQualifierValuePK.termId = :termId"), @NamedQuery(name = "DbxrefQualifierValueEntity.findByRank", query = "SELECT d FROM DbxrefQualifierValueEntity d WHERE d.dbxrefQualifierValuePK.rank = :rank"), @NamedQuery(name = "DbxrefQualifierValueEntity.findByValue", query = "SELECT d FROM DbxrefQualifierValueEntity d WHERE d.value = :value")})
public class DbxrefQualifierValueEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    protected DbxrefQualifierValueEntityPK dbxrefQualifierValuePK;
    @Column(name = "VALUE")
    private String value;
    @JoinColumn(name = "DBXREF_ID", referencedColumnName = "DBXREF_ID", insertable = false, updatable = false)
    @ManyToOne
    private DbxrefEntity dbxref;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID", insertable = false, updatable = false)
    @ManyToOne
    private TermEntity term;

    public DbxrefQualifierValueEntity() {
    }

    public DbxrefQualifierValueEntity(DbxrefQualifierValueEntityPK dbxrefQualifierValuePK) {
        this.dbxrefQualifierValuePK = dbxrefQualifierValuePK;
    }

    public DbxrefQualifierValueEntity(int dbxrefId, int termId, int rank) {
        this.dbxrefQualifierValuePK = new DbxrefQualifierValueEntityPK(dbxrefId, termId, rank);
    }

    public DbxrefQualifierValueEntityPK getDbxrefQualifierValuePK() {
        return dbxrefQualifierValuePK;
    }

    public void setDbxrefQualifierValuePK(DbxrefQualifierValueEntityPK dbxrefQualifierValuePK) {
        this.dbxrefQualifierValuePK = dbxrefQualifierValuePK;
    }

    public String getValue() {
        return value;
    }

    public void setValue(String value) {
        this.value = value;
    }

    public DbxrefEntity getDbxref() {
        return dbxref;
    }

    public void setDbxref(DbxrefEntity dbxref) {
        this.dbxref = dbxref;
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
        hash += (dbxrefQualifierValuePK != null ? dbxrefQualifierValuePK.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof DbxrefQualifierValueEntity)) {
            return false;
        }
        DbxrefQualifierValueEntity other = (DbxrefQualifierValueEntity) object;
        if ((this.dbxrefQualifierValuePK == null && other.dbxrefQualifierValuePK != null) || (this.dbxrefQualifierValuePK != null && !this.dbxrefQualifierValuePK.equals(other.dbxrefQualifierValuePK))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.DbxrefQualifierValueEntity[dbxrefQualifierValuePK=" + dbxrefQualifierValuePK + "]";
    }

}
