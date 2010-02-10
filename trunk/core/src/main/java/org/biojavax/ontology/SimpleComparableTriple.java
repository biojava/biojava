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

import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;

import org.biojava.bio.Annotation;
import org.biojava.ontology.Ontology;
import org.biojava.ontology.Term;
import org.biojava.ontology.Triple;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.RichAnnotation;

/**
 * Basic comparable triple, BioSQL style.
 * @author Richard Holland
 * @since 1.5
 */
public class SimpleComparableTriple extends AbstractChangeable implements ComparableTriple {
    
    private ComparableOntology ontology;
    private ComparableTerm object;
    private ComparableTerm subject;
    private ComparableTerm predicate;
    private Set descriptors = new TreeSet();
    
    /**
     * Creates a new instance of SimpleComparableTriple. All parameters are 
     * required and immutable.
     * @param ontology the ontology of the triple.
     * @param subject the subject of the triple.
     * @param object the object of the triple.
     * @param predicate the predicate of the triple.
     */
    SimpleComparableTriple(ComparableOntology ontology, ComparableTerm subject, ComparableTerm object, ComparableTerm predicate) {
        if (ontology == null) throw new IllegalArgumentException("Ontology must not be null");
        if (subject == null) throw new IllegalArgumentException("Subject must not be null");
        if (object == null) throw new IllegalArgumentException("Object must not be null");
        if (predicate == null) throw new IllegalArgumentException("Predicate must not be null");
        this.ontology = ontology;
        this.subject = subject;
        this.object = object;
        this.predicate = predicate;
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleComparableTriple() {}
    
    /**
     * {@inheritDoc}
     * Triples are sorted in order of ontology, subject, object, and finally
     * predicate.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        Triple them = (Triple)o;
        // Hibernate comparison - we haven't been populated yet
        if (this.ontology==null) return -1;
        // Normal comparison
        if (!this.ontology.equals(them.getOntology())) return this.ontology.compareTo((ComparableOntology)them.getOntology());
        if (!this.subject.equals(them.getSubject())) return this.subject.compareTo((ComparableTerm)them.getSubject());
        if (!this.object.equals(them.getObject())) return this.object.compareTo((ComparableTerm)them.getObject());
        return this.predicate.compareTo((ComparableTerm)them.getPredicate());
    }
    
    /**
     * {@inheritDoc}
     * Triples are equal only if they are from the same ontology and share the
     * same subject, object and predicate.
     */
    public boolean equals(Object o) {
        if(this == o) return true;
        if (o==null || !(o instanceof Triple)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.ontology==null) return false;
        // Normal comparison
        Triple them = (Triple)o;
        return (this.ontology.equals(them.getOntology()) &&
                this.subject.equals(them.getSubject()) &&
                this.object.equals(them.getObject()) &&
                this.predicate.equals(them.getPredicate()));
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int code = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.ontology==null) return code;
        // Normal comparison
        code = 37*code + this.ontology.hashCode();
        code = 37*code + this.subject.hashCode();
        code = 37*code + this.object.hashCode();
        code = 37*code + this.predicate.hashCode();
        return code;
    }
    
    /**
     * {@inheritDoc}
     * Returns the output of toSring()
     */
    public String getName() { return this.toString(); }
    
    /**
     * {@inheritDoc}
     */
    public Term getSubject() { return this.subject; }
    
    // Hibernate requirement - not for public use.
    void setSubject(ComparableTerm subject) { this.subject = subject; }
    
    /**
     * {@inheritDoc}
     */
    public Term getObject() { return this.object; }
    
    // Hibernate requirement - not for public use.
    void setObject(ComparableTerm object) { this.object = object; }
    
    /**
     * {@inheritDoc}
     */
    public Term getPredicate() { return this.predicate; }
    
    // Hibernate requirement - not for public use.
    void setPredicate(ComparableTerm predicate) { this.predicate = predicate; }
    
    /**
     * {@inheritDoc}
     */
    public void addDescriptor(ComparableTerm desc) throws IllegalArgumentException,ChangeVetoException {
        if (desc==null) throw new IllegalArgumentException("Cannot have null descriptor");
        if(!this.hasListeners(ComparableTriple.DESCRIPTOR)) {
            this.descriptors.add(desc);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    ComparableTriple.DESCRIPTOR,
                    desc,
                    null
                    );
            ChangeSupport cs = this.getChangeSupport(ComparableTriple.DESCRIPTOR);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.descriptors.add(desc);
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public boolean removeDescriptor(ComparableTerm desc) throws IllegalArgumentException,ChangeVetoException {
        if (desc==null) throw new IllegalArgumentException("Cannot have null descriptor");
        boolean result;
        if(!this.hasListeners(ComparableTriple.DESCRIPTOR)) {
            result = this.descriptors.remove(desc);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    ComparableTriple.DESCRIPTOR,
                    null,
                    desc
                    );
            ChangeSupport cs = this.getChangeSupport(ComparableTriple.DESCRIPTOR);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                result = this.descriptors.remove(desc);
                cs.firePostChangeEvent(ce);
            }
        }
        return result;
    }
    
    /**
     * {@inheritDoc}
     * <b>Warning</b> this method gives access to the original 
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public Set getDescriptors() { return this.descriptors; } // originals for Hibernate
    
    /**
     * {@inheritDoc}
     * <b>Warning</b> this method gives access to the original 
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public void setDescriptors(Set descriptors) throws ChangeVetoException {
        this.descriptors = descriptors;  // originals for Hibernate
    }
    
    /**
     * {@inheritDoc}
     * NOT IMPLEMENTED
     */
    public void removeSynonym(Object synonym) {
        throw new UnsupportedOperationException("BioJavaX does not know about triple synonyms.");
    }
    
    /**
     * {@inheritDoc}
     * NOT IMPLEMENTED
     */
    public void addSynonym(Object synonym) {
        throw new UnsupportedOperationException("BioJavaX does not know about triple synonyms.");
    }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS AN EMPTY LIST
     */
    public Object[] getSynonyms() { return Collections.EMPTY_LIST.toArray(); }
    
    /**
     * {@inheritDoc}
     */
    public Ontology getOntology() { return this.ontology; }
    
    // Hibernate requirement - not for public use.
    void setOntology(ComparableOntology ontology) { this.ontology = ontology; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY STRING
     */
    public String getDescription() { return ""; }
    
    /**
     * {@inheritDoc}
     * does not do anything
     */
    public void setDescription(String desc) { }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY ANNOTATION
     */
    public Annotation getAnnotation() { return RichAnnotation.EMPTY_ANNOTATION; }
    
    /**
     * {@inheritDoc}
     * Form: "ontology:predicate(subject,object)"
     */
    public String toString() {
        return this.ontology+":"+this.predicate+"("+this.subject+","+this.object+")";
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
