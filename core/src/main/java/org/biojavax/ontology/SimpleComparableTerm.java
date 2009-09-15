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

package org.biojavax.ontology;

import java.util.Arrays;
import java.util.Set;
import java.util.TreeSet;

import org.biojava.bio.Annotation;
import org.biojava.ontology.Ontology;
import org.biojava.ontology.Term;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.RankedCrossRef;
import org.biojavax.RichAnnotation;

/**
 * A Term object that can be compared and thus sorted.
 * @author Richard Holland
 * @since 1.5
 */
public class SimpleComparableTerm extends AbstractChangeable implements ComparableTerm {
    
    private String name;
    private String description;
    private ComparableOntology ontology;
    private String identifier;
    private Boolean obsolete;
    private Set synonyms = new TreeSet();
    private Set rankedcrossrefs = new TreeSet();
    
    /**
     * Creates a new instance of SimpleComparableTerm with synonyms.
     * @param ontology The ontology to put the term in. Must not be null.
     * @param name the name of the term. Must not be null.
     * @param synonyms a set of synonyms for the term. Can be null.
     */
    SimpleComparableTerm(ComparableOntology ontology, String name, Object[] synonyms) {
        if (name == null || name.equals("")) throw new IllegalArgumentException("Name must not be null or empty");
        if (ontology == null) throw new IllegalArgumentException("Ontology must not be null");        
        this.name = name;
        this.description = null;
        this.ontology = ontology;
        this.identifier = null;
        this.obsolete = Boolean.FALSE;
        if (synonyms!=null 
                && synonyms.length != 0) this.synonyms.addAll(Arrays.asList(synonyms));
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleComparableTerm() {}
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int value = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.ontology==null) return value;
        // Normal comparison
        value = 37*value + this.name.hashCode();
        value = 37*value + this.ontology.hashCode();
        return value;
    }
    
    /**
     * {@inheritDoc}
     * Two terms are equal if they are in the same ontology and
     * share the same name.
     */
    public boolean equals(Object obj) {
        if (obj == this) return true;
        if (!(obj instanceof Term)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.ontology==null) return false;
        // Normal comparison
        Term that = (Term) obj;
        return this.ontology.equals(that.getOntology()) &&
                this.name.equals(that.getName());
    }
    
    /**
     * {@inheritDoc}
     * Terms are sorted by ontology first, then name.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        // Hibernate comparison - we haven't been populated yet
        if (this.ontology==null) return -1;
        // Normal comparison
        Term them = (Term)o;
        if (!this.ontology.equals(them.getOntology())) return this.ontology.compareTo(them.getOntology());
        return this.name.compareTo(them.getName());
    }
    
    /**
     * {@inheritDoc}
     * Synonyms are stored in the database as the results of a toString() operation
     * on each synonym object. This doesn't happen until it reaches the database
     * though, so if you are not using a database, don't worry about it.
     */
    public void addSynonym(Object synonym) { this.synonyms.add(synonym); }
    
    /**
     * {@inheritDoc}
     */
    public void removeSynonym(Object synonym) { this.synonyms.remove(synonym); }
    
    /**
     * {@inheritDoc}
     */
    public Object[] getSynonyms() { return this.synonyms.toArray(); }
    
    // Hibernate requirement - not for public use.
    Set getSynonymSet() { return this.synonyms; }
    
    // Hibernate requirement - not for public use.
    void setSynonymSet(Set synonyms) { this.synonyms = synonyms; }
    
    /**
     * {@inheritDoc}
     * <b>Warning</b> this method gives access to the original 
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public Set getRankedCrossRefs() { return this.rankedcrossrefs; } // original for Hibernate
    
    /**
     * {@inheritDoc}
     * <b>Warning</b> this method gives access to the original 
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public void setRankedCrossRefs(Set rankedcrossrefs) throws ChangeVetoException {
        this.rankedcrossrefs = rankedcrossrefs; // original for Hibernate
    }
    
    /**
     * {@inheritDoc}
     */
    public void addRankedCrossRef(RankedCrossRef crossref) throws ChangeVetoException {
        if (crossref==null) throw new IllegalArgumentException("Crossref cannot be null");
        if(!this.hasListeners(ComparableTerm.RANKEDCROSSREF)) {
            this.rankedcrossrefs.add(crossref);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    ComparableTerm.RANKEDCROSSREF,
                    crossref,
                    null
                    );
            ChangeSupport cs = this.getChangeSupport(ComparableTerm.RANKEDCROSSREF);
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
        if(!this.hasListeners(ComparableTerm.RANKEDCROSSREF)) {
            this.rankedcrossrefs.remove(crossref);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    ComparableTerm.RANKEDCROSSREF,
                    null,
                    crossref
                    );
            ChangeSupport cs = this.getChangeSupport(ComparableTerm.RANKEDCROSSREF);
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
    public String getName() { return this.name; }
    
    // Hibernate requirement - not for public use.
    void setName(String name) { this.name = name; }
    
    /**
     * {@inheritDoc}
     */
    public String getDescription() { return this.description; }
    
    /**
     * {@inheritDoc}
     */
    public void setDescription(String description) throws ChangeVetoException {
        if(!this.hasListeners(ComparableTerm.DESCRIPTION)) {
            this.description = description;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    ComparableTerm.DESCRIPTION,
                    description,
                    this.description
                    );
            ChangeSupport cs = this.getChangeSupport(ComparableTerm.DESCRIPTION);
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
    public Ontology getOntology() { return this.ontology; }
    
    // Hibernate requirement - not for public use.
    void setOntology(ComparableOntology ontology) { this.ontology = ontology; }
    
    /**
     * {@inheritDoc}
     * Form: "ontology:name [obsolete]" where [obsolete] is optional
     */
    public String toString() { 
        boolean isobs = (this.obsolete!=null && this.obsolete.booleanValue());
        return this.ontology+":"+this.name+(isobs?" [obsolete]":""); 
    }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS AN EMPTY ANNOTATION OBJECT
     */
    public Annotation getAnnotation() { return RichAnnotation.EMPTY_ANNOTATION; }
    
    /**
     * {@inheritDoc}
     */
    public String getIdentifier() { return this.identifier; }
    
    /**
     * {@inheritDoc}
     */
    public void setIdentifier(String identifier) throws ChangeVetoException {
        if(!this.hasListeners(ComparableTerm.IDENTIFIER)) {
            this.identifier = identifier;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    ComparableTerm.IDENTIFIER,
                    identifier,
                    this.identifier
                    );
            ChangeSupport cs = this.getChangeSupport(ComparableTerm.IDENTIFIER);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.identifier = identifier;
                cs.firePostChangeEvent(ce);
            }
        }
    }

    // Hibernate requirement - not for public use.
    String getObsoleteChar() {
        return (this.getObsolete()!=null && this.getObsolete().equals(Boolean.TRUE))?"X":null;
    }

    // Hibernate requirement - not for public use.
    void setObsoleteChar(String obsolete) throws ChangeVetoException {
        this.setObsolete(Boolean.valueOf(obsolete!=null && obsolete.equals("X")));
    }
        
    /**
     * {@inheritDoc}
     */
    public Boolean getObsolete() { return this.obsolete; }
    
    /**
     * {@inheritDoc}
     */
    public void setObsolete(Boolean obsolete) throws ChangeVetoException {
        if(!this.hasListeners(ComparableTerm.OBSOLETE)) {
            this.obsolete = obsolete;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    ComparableTerm.OBSOLETE,
                    obsolete,
                    this.obsolete
                    );
            ChangeSupport cs = this.getChangeSupport(ComparableTerm.OBSOLETE);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.obsolete = obsolete;
                cs.firePostChangeEvent(ce);
            }
        }
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