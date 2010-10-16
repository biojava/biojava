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
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.RichObjectFactory;

/**
 * Reference implementation of NCBITaxon.
 * @author Richard Holland
 * @author Mark Schreiber
 * @author David Scott
 * @since 1.5
 */
public class SimpleNCBITaxon extends AbstractChangeable implements NCBITaxon {
    
    protected Set names = new TreeSet();
    private Map namesMap = new TreeMap();
    private Integer parent;
    private int NCBITaxID;
    private String nodeRank;
    private Integer geneticCode;
    private Integer mitoGeneticCode;
    private Integer leftValue;
    private Integer rightValue;
	private boolean isTaxonHidden = false;
	 
	private final static String TAXONISHIDDEN = "X";
	private final static int ROOTNCBIID = 1;
   
    /**
     * Creates a new instance of SimpleNCBITaxon based on the given taxon ID.
     * @param NCBITaxID the underlying taxon ID from NCBI.
     */
    public SimpleNCBITaxon(int NCBITaxID) {
        this.NCBITaxID = NCBITaxID;
    }
    
    /**
     * Creates a new instance of SimpleNCBITaxon based on the given taxon ID.
     * It may not be null, else you'll get exceptions.
     * @param NCBITaxID the underlying taxon ID from NCBI.
     */
    public SimpleNCBITaxon(Integer NCBITaxID) {
        this.NCBITaxID = NCBITaxID.intValue();
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleNCBITaxon() {}
    
    /**
     * {@inheritDoc}
     * NCBITaxon objects are compared only by their NCBITaxID fields.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        NCBITaxon them = (NCBITaxon)o;
        return this.NCBITaxID-them.getNCBITaxID();
    }
    
    /**
     * {@inheritDoc}
     * NCBITaxon objects are equal if their NCBITaxID fields match.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj==null || !(obj instanceof NCBITaxon)) return false;
        NCBITaxon them = (NCBITaxon)obj;
        return this.NCBITaxID==them.getNCBITaxID();
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        return this.NCBITaxID;
    }
    
    /**
     * {@inheritDoc}
     */
    public Set getNameClasses() { return this.namesMap.keySet(); }
    
    /**
     * {@inheritDoc}
     */
    public Set getNames(String nameClass) throws IllegalArgumentException {
        if (nameClass==null) throw new IllegalArgumentException("Name class cannot be null");
        Set items = (Set)this.namesMap.get(nameClass);
        Set n = new TreeSet();
        if (items!=null) for (Iterator j = items.iterator(); j.hasNext(); ) {
            SimpleNCBITaxonName name = (SimpleNCBITaxonName)j.next();
            n.add(name.getName());
        }
        return n;
    }
    
    // Hibernate requirement - not for public use.
    protected Set getNameSet() { return this.names; } // original for Hibernate
    
    // Hibernate requirement - not for public use.
    void setNameSet(Set names) {
        this.names = names; // original for Hibernate
        // convert set to map
        this.namesMap.clear();
        for (Iterator i = names.iterator(); i.hasNext(); ) {
            SimpleNCBITaxonName n = (SimpleNCBITaxonName)i.next();
            try {
                this.addName(n.getNameClass(), n.getName());
            } catch (ChangeVetoException e) {
                throw new RuntimeException("Database contents don't add up",e);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void addName(String nameClass, String name) throws IllegalArgumentException,ChangeVetoException {
        if (name==null) throw new IllegalArgumentException("Name cannot be null");
        if (nameClass==null) throw new IllegalArgumentException("Name class cannot be null");
        SimpleNCBITaxonName n = new SimpleNCBITaxonName(nameClass, name);
        if(!this.hasListeners(NCBITaxon.NAMES)) {
            if (!this.namesMap.containsKey(nameClass)) this.namesMap.put(nameClass,new TreeSet());
            ((Set)this.namesMap.get(nameClass)).add(n);
            this.names.add(n);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    NCBITaxon.NAMES,
                    name,
                    ((Set)this.namesMap.get(nameClass)).contains(n)?name:null
                    );
            ChangeSupport cs = this.getChangeSupport(NCBITaxon.NAMES);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                if (!this.namesMap.containsKey(nameClass)) this.namesMap.put(nameClass,new TreeSet());
                ((Set)this.namesMap.get(nameClass)).add(n);
                this.names.add(n);
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public boolean removeName(String nameClass, String name) throws IllegalArgumentException,ChangeVetoException {
        if (name==null) throw new IllegalArgumentException("Name cannot be null");
        if (nameClass==null) throw new IllegalArgumentException("Name class cannot be null");
        SimpleNCBITaxonName n = new SimpleNCBITaxonName(nameClass, name);
        if (!this.namesMap.containsKey(nameClass)) return false;
        boolean results;
        if(!this.hasListeners(NCBITaxon.NAMES)) {
            results = ((Set)this.namesMap.get(nameClass)).remove(n);
            this.names.remove(n);
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    NCBITaxon.NAMES,
                    null,
                    name
                    );
            ChangeSupport cs = this.getChangeSupport(NCBITaxon.NAMES);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                results = ((Set)this.namesMap.get(nameClass)).remove(n);
                this.names.remove(n);
                cs.firePostChangeEvent(ce);
            }
        }
        return results;
    }
    
    /**
     * {@inheritDoc}
     */
    public boolean containsName(String nameClass, String name) throws IllegalArgumentException {
        if (name==null) throw new IllegalArgumentException("Name cannot be null");
        if (nameClass==null) throw new IllegalArgumentException("Name class cannot be null");
        if (!this.namesMap.containsKey(nameClass)) return false;
        SimpleNCBITaxonName n = new SimpleNCBITaxonName(nameClass, name);
        return ((Set)this.namesMap.get(nameClass)).contains(n);
    }
    
    protected final Map getNamesMap() {
    	return namesMap;
    }
    
    /**
     * {@inheritDoc}
     */
    public Integer getParentNCBITaxID() { return this.parent; }
    
    /**
     * {@inheritDoc}
     */
    public void setParentNCBITaxID(Integer parent) throws ChangeVetoException {
        if(!this.hasListeners(NCBITaxon.PARENT)) {
            this.parent = parent;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    NCBITaxon.PARENT,
                    parent,
                    this.parent
                    );
            ChangeSupport cs = this.getChangeSupport(NCBITaxon.PARENT);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.parent = parent;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public int getNCBITaxID() { return this.NCBITaxID; }
    
    // Hibernate requirement - not for public use.
    void setNCBITaxID(int NCBITaxID) { this.NCBITaxID = NCBITaxID; }
    
    /**
     * {@inheritDoc}
     */
    public String getNodeRank() { return this.nodeRank; }
    
    /**
     * Setter for property nodeRank.
     * @param nodeRank New value of property nodeRank.
     * @throws org.biojava.utils.ChangeVetoException in case of objections.
     */
    public void setNodeRank(String nodeRank) throws ChangeVetoException {
        if(!this.hasListeners(NCBITaxon.NODERANK)) {
            this.nodeRank = nodeRank;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    NCBITaxon.NODERANK,
                    nodeRank,
                    this.nodeRank
                    );
            ChangeSupport cs = this.getChangeSupport(NCBITaxon.NODERANK);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.nodeRank = nodeRank;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public Integer getGeneticCode() { return this.geneticCode; }
    
    /**
     * {@inheritDoc}
     */
    public void setGeneticCode(Integer geneticCode) throws ChangeVetoException {
        if(!this.hasListeners(NCBITaxon.GENETICCODE)) {
            this.geneticCode = geneticCode;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    NCBITaxon.GENETICCODE,
                    nodeRank,
                    this.nodeRank
                    );
            ChangeSupport cs = this.getChangeSupport(NCBITaxon.GENETICCODE);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.geneticCode = geneticCode;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * Getter for property mitoGeneticCode. Returns Persistent.NULL_INTEGER if null.
     * @return Value of property mitoGeneticCode.
     */
    public Integer getMitoGeneticCode() { return this.mitoGeneticCode; }
    
    /**
     * {@inheritDoc}
     */
    public void setMitoGeneticCode(Integer mitoGeneticCode) throws ChangeVetoException {
        if(!this.hasListeners(NCBITaxon.MITOGENETICCODE)) {
            this.mitoGeneticCode = mitoGeneticCode;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    NCBITaxon.MITOGENETICCODE,
                    mitoGeneticCode,
                    this.mitoGeneticCode
                    );
            ChangeSupport cs = this.getChangeSupport(NCBITaxon.MITOGENETICCODE);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.mitoGeneticCode = mitoGeneticCode;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public Integer getLeftValue() { return this.leftValue; }
    
    /**
     * {@inheritDoc}
     */
    public void setLeftValue(Integer leftValue) throws ChangeVetoException {
        if(!this.hasListeners(NCBITaxon.LEFTVALUE)) {
            this.leftValue = leftValue;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    NCBITaxon.LEFTVALUE,
                    leftValue,
                    this.leftValue
                    );
            ChangeSupport cs = this.getChangeSupport(NCBITaxon.LEFTVALUE);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.leftValue = leftValue;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public Integer getRightValue() { return this.rightValue; }
    
    /**
     * {@inheritDoc}
     */
    public void setRightValue(Integer rightValue) throws ChangeVetoException {
        if(!this.hasListeners(NCBITaxon.RIGHTVALUE)) {
            this.rightValue = rightValue;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    NCBITaxon.RIGHTVALUE,
                    rightValue,
                    this.rightValue
                    );
            ChangeSupport cs = this.getChangeSupport(NCBITaxon.RIGHTVALUE);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.rightValue = rightValue;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * Returns the name of this taxon entry in the form:
     *   scientific (common)
     * or if there is no common name:
     *   scientific
     * or if there are no scientific names at all, the empty string.
     * @return the display name as described above.
     */
    public String getDisplayName() {
        StringBuffer sb = new StringBuffer();
        Set sciNames = this.getNames(NCBITaxon.SCIENTIFIC);
        Set comNames = this.getNames(NCBITaxon.COMMON);
        if (sciNames.size()>0) {
            sb.append((String)sciNames.iterator().next());
            if (comNames.size()>0) {
                sb.append(" (");
                sb.append((String)comNames.iterator().next());
                sb.append(")");
            }
        }
        return sb.toString();
    }
    
    /**
     * Returns the taxonomy hierarchy of this taxon entry in the form:
     *   most specific; less specific; ...; least specific.
     * It follows the chain up the tree as far as it can, and will stop
     * at the first one it comes to that returns null for getParentNCBITaxID().
     * If this taxon entry has no scientific name, you will get the string ".".
     * @return the display name as described above.
     */
    /*
     * (non-Javadoc)
     * @see com.bioperception.bio.taxa.NCBITaxon2#isTaxonHidden()
     */
    public final boolean isTaxonHidden() {
    	return isTaxonHidden;
    }
    
    /*
     * (non-Javadoc)
     * @see com.bioperception.bio.taxa.NCBITaxon2#setTaxonHidden(boolean)
     */
    public final void setTaxonHidden(final boolean isTaxonHidden) throws ChangeVetoException {
        if(!this.hasListeners(HIDDEN)) {
            this.isTaxonHidden = isTaxonHidden;
        } else {
            final ChangeEvent ce = new ChangeEvent( this, HIDDEN, new Boolean(isTaxonHidden), new Boolean(this.isTaxonHidden));
            final ChangeSupport cs = this.getChangeSupport(HIDDEN);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.isTaxonHidden = isTaxonHidden;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    // Hibernate requirement - not for public use.
    String getTaxonHiddenChar() {
        return isTaxonHidden()?TAXONISHIDDEN:null;
    }

    // Hibernate requirement - not for public use.
    void setTaxonHiddenChar(final String isHiddenChar) throws ChangeVetoException {
        setTaxonHidden(isHiddenChar!=null || (isHiddenChar!=null && isHiddenChar.length() > 0));// any character will set
    }

    /**
     * Returns the taxonomy hierarchy of this taxon entry in the form:
     *   most specific; less specific; ...; least specific.
     * It follows the chain up the tree as far as it can, and will stop
     * at the first one it comes to that returns null for getParentNCBITaxID().
     * If this taxon entry has no scientific name, you will get the string ".".
     * @return the display name as described above.
     */
    public String getNameHierarchy() {
    	StringBuffer sb = new StringBuffer();
    	boolean first = true;
    	Integer parent = this.getParentNCBITaxID();
    	while (parent!=null) {
    		NCBITaxon t = (NCBITaxon)RichObjectFactory.getObject(SimpleNCBITaxon.class, new Object[]{parent});
    		Set sciNames = t.getNames(NCBITaxon.SCIENTIFIC);
//    		System.out.println("SimpleNCBITaxon2.getNameHierarchy-t:"+t+", t.isTaxonHidden? "+t.isTaxonHidden()+", isRoot? "+isRoot(t)+", sciNames:"+sciNames);
    		if (sciNames.size()>0) {
    			if (!t.isTaxonHidden() && !isRoot(t)) {// root is NOT hidden - but also not displayed
    				if (!first) sb.insert(0,"; ");
    				else first = false;
    				sb.insert(0,(String)sciNames.iterator().next());
    			}
    			// Don't get into endless loop if child's parent is itself.
    			if (t.getParentNCBITaxID().equals(new Integer(t.getNCBITaxID()))) parent = null;
    			else parent = t.getParentNCBITaxID();
    		}
    		// Also don't go up past a parent that doesn't have a name.
    		else parent = null;
    	}
    	sb.append(".");
    	return sb.toString();
    }
    
    
    private final static boolean isRoot(final NCBITaxon theTaxon) {
    	return isRoot(theTaxon.getNCBITaxID());
    }
    
    private final static boolean isRoot(final int theNCBIId) {
    	return theNCBIId==ROOTNCBIID;
    }
    
    /**
     * {@inheritDoc}
     * Form: "taxid:[name,name...]"
     */
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("taxid:");
        sb.append(this.NCBITaxID);
        sb.append("[");
        for (Iterator i = this.getNameSet().iterator(); i.hasNext(); ) {
            SimpleNCBITaxonName n = (SimpleNCBITaxonName)i.next();
            sb.append(n);
            if (i.hasNext()) sb.append(",");
        }
        sb.append("]");
        return sb.toString();
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

