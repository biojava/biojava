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

package org.biojavax.bio.taxa;

/**
 * Represents a name class plus name combination for an NCBITaxon object.
 * @author Richard Holland
 * @since 1.5
 */
public class SimpleNCBITaxonName implements Comparable {
    
    private String nameClass;
    private String name;

    // Hibernate requirement - not for public use.
    protected SimpleNCBITaxonName() {}
    
    /**
     * Creates a new taxon name based on the given class and name, both of
     * which cannot be null.
     * @param nameClass the class of the new name. Use one of the constants from
     * {@link org.biojavax.bio.taxa.NCBITaxon} (for example {@link org.biojavax.bio.taxa.NCBITaxon#SCIENTIFIC}).
     * @param name the name itself
     */
    public SimpleNCBITaxonName(String nameClass, String name) {
        if (nameClass==null) throw new IllegalArgumentException("Name class cannot be null");
        if (name==null) throw new IllegalArgumentException("Name cannot be null");
        if (name.indexOf('\n') >= 0) throw new IllegalArgumentException("NCBI taxonomy names cannot embed new lines - at:"+name.indexOf('\n')+", in name: <"+name+">");
        this.nameClass = nameClass;
        this.name = name; 
    }
    
    /**
     * Changes the class of this name.
     * @param nameClass the new class for this name.
     */
    public void setNameClass(String nameClass) { 
        if (nameClass==null) throw new IllegalArgumentException("Name class cannot be null");
        this.nameClass = nameClass; 
    }
    
    /**
     * Returns the class of this name.
     * @return the class of this name.
     */
    public String getNameClass() { return this.nameClass; }
    
    /**
     * Changes the name.
     * @param name the new name.
     */
    public void setName(String name) {       
        if (name==null) throw new IllegalArgumentException("Name cannot be null");
        this.name = name; 
    }
    
    /**
     * Returns this name.
     * @return this name.
     */
    public String getName() { return this.name; }
    
    /**
     * {@inheritDoc}
     * Two taxon names are equal if their name and class match.
     */
    public boolean equals(Object o) {
        if (o==this) return true;
        if (!(o instanceof SimpleNCBITaxonName)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.nameClass==null) return false;
        // Normal comparison
        SimpleNCBITaxonName them = (SimpleNCBITaxonName) o;
        return them.getNameClass().equals(this.nameClass) &&
                them.getName().equals(this.name);
    }    
    
    /**
     * {@inheritDoc}
     * Taxon names are sorted by class first, then name.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        // Hibernate comparison - we haven't been populated yet
        if (this.nameClass==null) return -1;
        // Normal comparison
        SimpleNCBITaxonName them = (SimpleNCBITaxonName)o;
        if (!them.getNameClass().equals(this.nameClass)) return this.nameClass.compareTo(them.getNameClass());
        return this.name.compareTo(them.getName());
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int code = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.nameClass==null) return code;
        // Normal comparison
        code = 31*code + this.name.hashCode();
        code = 31*code + this.nameClass.hashCode();
        return code;
    }
    
    /**
     * {@inheritDoc}
     * Form: "class:name"
     */
    public String toString() {
        return this.nameClass+":"+this.name;
    }
}
