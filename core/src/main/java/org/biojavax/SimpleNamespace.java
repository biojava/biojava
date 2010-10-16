/*
 * BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence. This should
 * be distributed with the code. If you do not have a copy,
 * see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors. These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 * http://www.biojava.org/
 *
 */

package org.biojavax;

import java.net.URI;
import java.net.URISyntaxException;

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * A basic Namespace implemenation.
 * @author Richard Holland
 * @author Mark Schreiber
 * @since 1.5
 */

public class SimpleNamespace extends AbstractChangeable implements Namespace {
    
    private String name;
    private String acronym;
    private String authority;
    private String description;
    private URI URI;
    
    /**
     * Creates a new instance of SimpleNamespace with the given name,
     * which cannot be null.
     * @param name the name of the namespace.
     */
    public SimpleNamespace(String name) {
        if (name==null) throw new IllegalArgumentException("Name cannot be null");
        this.name = name;
        this.acronym = null;
        this.authority = null;
        this.description = null;
        this.URI = null;
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleNamespace() {}
    
    /**
     * {@inheritDoc}
     */
    public void setAcronym(String acronym) throws ChangeVetoException {
        if(!this.hasListeners(Namespace.ACRONYM)) {
            this.acronym = acronym;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    Namespace.ACRONYM,
                    acronym,
                    this.acronym
                    );
            ChangeSupport cs = this.getChangeSupport(Namespace.ACRONYM);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.acronym = acronym;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void setAuthority(String authority) throws ChangeVetoException {
        if(!this.hasListeners(Namespace.AUTHORITY)) {
            this.authority = authority; } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    Namespace.AUTHORITY,
                    authority,
                    this.authority
                    );
            ChangeSupport cs = this.getChangeSupport(Namespace.AUTHORITY);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.authority = authority;
                cs.firePostChangeEvent(ce);
            }
            }
    }
    
    /**
     * {@inheritDoc}
     */
    public void setDescription(String description) throws ChangeVetoException {
        if(!this.hasListeners(Namespace.DESCRIPTION)) {
            this.description = description;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    Namespace.DESCRIPTION,
                    description,
                    this.description
                    );
            ChangeSupport cs = this.getChangeSupport(Namespace.DESCRIPTION);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.description = description;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    // Hibernate requirement - not for public use.
    // Converts a String object representing a URI into an actual URI object.
    void setURIString(String URI) throws ChangeVetoException, URISyntaxException {
        if (URI!=null) this.setURI(new URI(URI));
        else this.URI=null;
    }
    
    // Hibernate requirement - not for public use.
    // Converts a URI object into a String representation of that URI
    String getURIString() {
        if (this.URI==null) return null;
        else return this.URI.toASCIIString();
    }
    
    /**
     * {@inheritDoc}
     */
    public void setURI(URI URI) throws ChangeVetoException {
        if(!this.hasListeners(Namespace.URI)) {
            this.URI = URI;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    Namespace.URI,
                    URI,
                    this.URI
                    );
            ChangeSupport cs = this.getChangeSupport(Namespace.URI);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.URI = URI;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    // Hibernate requirement - not for public use.
    void setName(String name) { this.name = name; }
    
    /**
     * {@inheritDoc}
     */
    public String getAcronym() { return this.acronym; }
    
    /**
     * {@inheritDoc}
     */
    public String getAuthority() { return this.authority; }
    
    /**
     * {@inheritDoc}
     */
    public String getDescription() { return this.description; }
    
    /**
     * {@inheritDoc}
     */
    public String getName() { return this.name; }
    
    /**
     * {@inheritDoc}
     */
    public URI getURI() { return this.URI; }
    
    /**
     * {@inheritDoc}
     * Namespaces are compared only by name.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        // Hibernate comparison - we haven't been populated yet
        if (this.name==null) return -1;
        // Normal comparison
        Namespace them = (Namespace)o;
        return this.name.compareTo(them.getName());
    }
    
    /**
     * {@inheritDoc}
     * Namespaces are equal only by name.
     */
    public boolean equals(Object obj) {
        if(this == obj) return true;
        if (obj==null || !(obj instanceof Namespace)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.name==null) return false;
        // Normal comparison
        Namespace them = (Namespace)obj;
        return this.name.equals(them.getName());
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int hash = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.name==null) return hash;
        // Normal comparison
        return 31*hash + this.name.hashCode();
    }
    
    /**
     * {@inheritDoc}
     * Form: "name"
     */
    public String toString() { 
        return this.getName(); 
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