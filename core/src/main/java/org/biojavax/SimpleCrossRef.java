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

package org.biojavax;

import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.utils.Unchangeable;

/**
 * A basic CrossRef implementation.
 * @author Richard Holland
 * @author Mark Schreiber
 * @author George Waldon (made Unchangeable)
 * @since 1.5
 */
public class SimpleCrossRef extends Unchangeable implements CrossRef {
    
    private RichAnnotation notes = new SimpleRichAnnotation();
    private String accession;
    private String dbname;
    private int version;
    
    /**
     * Creates a new instance of SimpleCrossRef with the values to use for
     * the immutable database name, accession and version.
     * @param dbname the dbname for this crossref.
     * @param accession the accession for this crossref.
     * @param version the version for this crossref.
     */
    public SimpleCrossRef(String dbname, String accession, int version) {
        if (dbname==null) throw new IllegalArgumentException("Database name cannot be null");
        if (accession==null) throw new IllegalArgumentException("Accession cannot be null");
        this.accession = accession;
        this.dbname = dbname;
        this.version = version;
    }
    
    /**
     * Creates a new instance of SimpleCrossRef with the values to use for
     * the immutable database name, accession and version. Identical to other
     * dbname/accession/version constructor except the version is specified
     * as an Integer object rather than an int primitive. Will throw an
     * exception if version is null.
     * @param dbname the dbname for this crossref.
     * @param accession the accession for this crossref.
     * @param version the version for this crossref.
     */
    public SimpleCrossRef(String dbname, String accession, Integer version) {
        this(dbname,accession,version.intValue());
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleCrossRef() {}
    
    /**
     * {@inheritDoc}
     */
    public Annotation getAnnotation() { return getRichAnnotation(); }

    /**
     * {@inheritDoc}
     */
    public RichAnnotation getRichAnnotation() { return this.notes; }

    /**
     * {@inheritDoc}
     */
    public Set getNoteSet() { return this.notes.getNoteSet(); }
    
    /**
     * {@inheritDoc}
     */
    public void setNoteSet(Set notes) { this.notes.setNoteSet(notes); }
    
    /**
     * {@inheritDoc}
     */
    public String getAccession() { return this.accession; }
    
    // Hibernate requirement - not for public use.
    void setAccession(String accession) { this.accession = accession; }
    
    /**
     * {@inheritDoc}
     */
    public String getDbname() { return this.dbname; }
    
    // Hibernate requirement - not for public use.
    void setDbname(String dbname) { this.dbname = dbname; }
    
    /**
     * {@inheritDoc}
     */
    public int getVersion() { return this.version; }
    
    // Hibernate requirement - not for public use.
    void setVersion(int version) { this.version = version; }
    
    /**
     * {@inheritDoc}
     * Compares cross references first by database name, then by accession,
     * then by version.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        // Hibernate comparison - we haven't been populated yet
        if (this.dbname==null) return -1;
        // Normal comparison
        CrossRef them = (CrossRef)o;
        if (!this.dbname.equals(them.getDbname())) return this.dbname.compareTo(them.getDbname());
        if (!this.accession.equals(them.getAccession())) return this.accession.compareTo(them.getAccession());
        return this.version-them.getVersion();
    }
    
    /**
     * {@inheritDoc}
     * Equality is defined as having the same database name, accession and
     * version.
     */
    public boolean equals(Object obj) {
        if(this == obj) return true;
        if (obj==null || !(obj instanceof CrossRef)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.dbname==null) return false;
        // Normal comparison
        CrossRef them = (CrossRef)obj;
        return (this.dbname.equals(them.getDbname()) &&
                this.accession.equals(them.getAccession()) &&
                this.version==them.getVersion()
                );
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int code = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.dbname==null) return code;
        // Normal comparison
        code = 37*code + this.dbname.hashCode();
        code = 37*code + this.accession.hashCode();
        code = 37*code + this.version;
        return code;
    }
    
    /**
     * {@inheritDoc}
     * Form: "dbname:accession.version"
     */
    public String toString() {
        return this.getDbname()+":"+this.getAccession()+"."+this.getVersion();
    }
    
    // Hibernate requirement - not for public use.
    private Integer id;
    
    /**
     * Gets the Hibernate ID. Should be used with caution.
     * @return the Hibernate ID, if using Hibernate.
     */
    public Integer getId() { return this.id; }
    
    /**
     * Sets the Hibernate ID. Should be used with caution.
     * @param id the Hibernate ID, if using Hibernate.
     */
    public void setId(Integer id) { this.id = id;}
}

