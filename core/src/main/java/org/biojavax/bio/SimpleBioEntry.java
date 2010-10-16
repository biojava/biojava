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

package org.biojavax.bio;
import java.util.Set;
import java.util.TreeSet;

import org.biojava.bio.Annotation;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.Comment;
import org.biojavax.Namespace;
import org.biojavax.RankedCrossRef;
import org.biojavax.RankedDocRef;
import org.biojavax.RichAnnotation;
import org.biojavax.SimpleRichAnnotation;
import org.biojavax.bio.taxa.NCBITaxon;

/**
 * Reference implementation of a BioEntry object which has no features or sequence.
 * Equality is the combination of namespace, name, accession and version.
 * @author Richard Holland
 * @author Mark Schreiber
 * @author George Waldon
 * @since 1.5
 */
public class SimpleBioEntry extends AbstractChangeable implements BioEntry {
    
    private Set comments = new TreeSet();
    private Set rankedcrossrefs = new TreeSet();
    private Set rankeddocrefs = new TreeSet();
    private Set relationships = new TreeSet();
    private String description;
    private String division;
    private String identifier;
    private String name;
    private String accession;
    private int version;
    private NCBITaxon taxon;
    private Namespace ns;
    private RichAnnotation notes = new SimpleRichAnnotation();
    
    /**
     * Creates a new bioentry representing the sequence in the given namespace
     * with the given name, accession and version. These properties are all
     * immutable and non-nullable.
     * @param ns The namespace for this new bioentry (not null).
     * @param name The name for this new bioentry (not null).
     * @param accession The accession for this new bioentry (not null).
     * @param version The version for this new bioentry.
     */
    public SimpleBioEntry(Namespace ns, String name, String accession, int version) {
        if (name==null) throw new IllegalArgumentException("Name cannot be null");
        if (accession==null) throw new IllegalArgumentException("Accession cannot be null");
        if (ns==null) throw new IllegalArgumentException("Namespace cannot be null");
        this.description = null;
        this.division = null;
        this.identifier = null;
        this.name = name;
        this.accession = accession;
        this.version = version;
        this.taxon = null;
        this.ns = ns;
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleBioEntry() {} // protected so SimpleRichSequence can extend us
    
    /**
     * {@inheritDoc} 
     * <b>Warning</b> this method gives access to the original 
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public Set getRankedCrossRefs() { return this.rankedcrossrefs; } // original for Hibernate
    
    /**
     * {@inheritDoc}
     */
    public void setTaxon(NCBITaxon taxon) throws ChangeVetoException {
        if(!this.hasListeners(BioEntry.TAXON)) {
            this.taxon = taxon;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.TAXON,
                    taxon,
                    this.taxon
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.TAXON);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.taxon = taxon;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
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
     * <b>Warning</b> this method gives access to the original 
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public Set getNoteSet() { return this.notes.getNoteSet(); }
    
    /**
     * {@inheritDoc}
     * <b>Warning</b> this method gives access to the original 
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public void setNoteSet(Set notes) throws ChangeVetoException { this.notes.setNoteSet(notes); }
    
    /**
     * {@inheritDoc}
     * <b>Warning</b> this method gives access to the original 
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public Set getComments() { return this.comments; } // must be original for Hibernate
    
    /**
     * {@inheritDoc}
     * <b>Warning</b> this method gives access to the original 
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public Set getRankedDocRefs() { return this.rankeddocrefs; } // must be original for Hibernate
    
    /**
     * {@inheritDoc}
     * <b>Warning</b> this method gives access to the original 
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public Set getRelationships() { return this.relationships; }  // must be original for Hibernate
    
    /**
     * {@inheritDoc}
     */
    public void setIdentifier(String identifier) throws ChangeVetoException {
        if(!this.hasListeners(BioEntry.IDENTIFIER)) {
            this.identifier = identifier;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.IDENTIFIER,
                    identifier,
                    this.identifier
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.IDENTIFIER);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.identifier = identifier;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void setDivision(String division) throws ChangeVetoException {
        if(!this.hasListeners(BioEntry.DIVISION)) {
            this.division = division;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.DIVISION,
                    division,
                    this.division
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.DIVISION);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.division = division;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void setDescription(String description) throws ChangeVetoException {
        if(!this.hasListeners(BioEntry.DESCRIPTION)) {
            this.description = description;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.DESCRIPTION,
                    description,
                    this.description
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.DESCRIPTION);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.description = description;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public String getAccession() { return this.accession; }
    
    /**
     * {@inheritDoc}
     */
    public String getDescription() { return this.description; }
    
    /**
     * {@inheritDoc}
     */
    public String getDivision() { return this.division; }
    
    /**
     * {@inheritDoc}
     */
    public String getIdentifier() { return this.identifier; }
    
    /**
     * {@inheritDoc}
     */
    public String getName() { return this.name; }
    
    /**
     * {@inheritDoc}
     */
    public Namespace getNamespace() { return this.ns; }
    
    /**
     * {@inheritDoc}
     */
    public NCBITaxon getTaxon() { return this.taxon; }
    
    /**
     * {@inheritDoc}
     */
    public int getVersion() { return this.version; }
    
    /**
     * {@inheritDoc}
     * Two bioentries are equal if they share the same namespace, name,
     * accession and version.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj==null || !(obj instanceof BioEntry)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.ns==null) return false;
        // Normal comparison
            BioEntry them = (BioEntry)obj;
            return (this.ns.equals(them.getNamespace()) &&
                    this.name.equals(them.getName()) &&
                    this.accession.equals(them.getAccession()) &&
                    this.version==them.getVersion());
    }
    
    /**
     * {@inheritDoc}
     * Bioentries are ordered first by namespace, then name, accession, and
     * finally version.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        // Hibernate comparison - we haven't been populated yet
        if (this.ns==null) return -1;
        // Normal comparison
        BioEntry them = (BioEntry)o;
        if (!this.ns.equals(them.getNamespace())) return this.ns.compareTo(them.getNamespace());
        if (!this.name.equals(them.getName())) return this.name.compareTo(them.getName());
        if (!this.accession.equals(them.getAccession())) return this.accession.compareTo(them.getAccession());
        return this.version-them.getVersion();
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int code = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.ns==null) return code;
        // Normal comparison
        code = 37*code + this.ns.hashCode();
        code = 37*code + this.name.hashCode();
        code = 37*code + this.accession.hashCode();
        code = 37*code + this.version;
        return code;
    }
    
    /**
     * {@inheritDoc} 
     * Form: namespace:name/accession.version
     */
    public String toString() { 
        return this.getNamespace()+":"+this.getName()+"/"+this.getAccession()+"."+this.getVersion(); 
    }
        
    /**
     * {@inheritDoc}
     */
    public void addRankedCrossRef(RankedCrossRef crossref) throws ChangeVetoException {
        if (crossref==null) throw new IllegalArgumentException("Crossref cannot be null");
        if(!this.hasListeners(BioEntry.RANKEDCROSSREF)) {
            this.rankedcrossrefs.add(crossref);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.RANKEDCROSSREF,
                    crossref,
                    null
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.RANKEDCROSSREF);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.rankedcrossrefs.add(crossref);
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void removeRankedCrossRef(RankedCrossRef crossref) throws ChangeVetoException {
        if (crossref==null) throw new IllegalArgumentException("Crossref cannot be null");
        if(!this.hasListeners(BioEntry.RANKEDCROSSREF)) {
            this.rankedcrossrefs.remove(crossref);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.RANKEDCROSSREF,
                    null,
                    crossref
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.RANKEDCROSSREF);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.rankedcrossrefs.remove(crossref);
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void addRankedDocRef(RankedDocRef docref) throws ChangeVetoException {
        if (docref==null) throw new IllegalArgumentException("Docref cannot be null");
        if(!this.hasListeners(BioEntry.RANKEDDOCREF)) {
            this.rankeddocrefs.add(docref);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.RANKEDDOCREF,
                    docref,
                    null
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.RANKEDDOCREF);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.rankeddocrefs.add(docref);
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void removeRankedDocRef(RankedDocRef docref) throws ChangeVetoException {
        if (docref==null) throw new IllegalArgumentException("Docref cannot be null");
        if(!this.hasListeners(BioEntry.RANKEDDOCREF)) {
            this.rankeddocrefs.remove(docref);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.RANKEDDOCREF,
                    null,
                    docref
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.RANKEDDOCREF);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.rankeddocrefs.remove(docref);
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void addComment(Comment comment) throws ChangeVetoException {
        if (comment==null) throw new IllegalArgumentException("Comment cannot be null");
        if(!this.hasListeners(BioEntry.COMMENT)) {
            this.comments.add(comment);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.COMMENT,
                    comment,
                    null
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.COMMENT);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.comments.add(comment);
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void removeComment(Comment comment) throws ChangeVetoException {
        if (comment==null) throw new IllegalArgumentException("Comment cannot be null");
        if(!this.hasListeners(BioEntry.COMMENT)) {
            this.comments.remove(comment);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.COMMENT,
                    null,
                    comment
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.COMMENT);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.comments.remove(comment);
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void addRelationship(BioEntryRelationship relation) throws ChangeVetoException {
        if (relation==null) throw new IllegalArgumentException("Relationship cannot be null");
        if(!this.hasListeners(BioEntry.RELATIONS)) {
            this.relationships.add(relation);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.RELATIONS,
                    relation,
                    null
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.RELATIONS);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.relationships.add(relation);
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void removeRelationship(BioEntryRelationship relation) throws ChangeVetoException {
        if (relation==null) throw new IllegalArgumentException("Relationship cannot be null");
        if(!this.hasListeners(BioEntry.RELATIONS)) {
            this.relationships.remove(relation);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntry.RELATIONS,
                    null,
                    relation
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntry.RELATIONS);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.relationships.remove(relation);
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    // Hibernate requirement - not for public use.
    void setRelationships(Set relationships) { this.relationships = relationships; }  // must be original for Hibernate
    
    // Hibernate requirement - not for public use.
    void setNamespace(Namespace ns) { this.ns = ns; }
    
    // Hibernate requirement - not for public use.
    void setName(String name) { this.name = name; }
    
    // Hibernate requirement - not for public use.
    void setAccession(String acc) { this.accession = acc; }
    
    // Hibernate requirement - not for public use.
    void setVersion(int v) { this.version = v; }
    
    // Hibernate requirement - not for public use.
    void setRankedDocRefs(Set docrefs) { this.rankeddocrefs = docrefs; } // must be original for Hibernate
    
    // Hibernate requirement - not for public use.
    void setComments(Set comments) { this.comments = comments; }  // must be original for Hibernate
    
    // Hibernate requirement - not for public use.
    public void setRankedCrossRefs(Set rankedcrossrefs) { this.rankedcrossrefs = rankedcrossrefs; } // original for Hibernate
    
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

