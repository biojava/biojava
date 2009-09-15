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
import java.util.Set;

import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * Represents an NCBI Taxon entry, a combination of the taxon and taxon_name
 * tables in BioSQL.
 * @author Mark Schreiber
 * @author Richard Holland
 * @since 1.5
 */
public interface NCBITaxon extends Comparable,Changeable {
    
    public static final ChangeType NAMES = new ChangeType(
            "This taxon's names have changed",
            "org.biojavax.bio.taxa.NCBITaxon",
            "NAMES"
            );
    public static final ChangeType PARENT = new ChangeType(
            "This taxon's parent has changed",
            "org.biojavax.bio.taxa.NCBITaxon",
            "PARENT"
            );
    public static final ChangeType NODERANK = new ChangeType(
            "This taxon's node rank has changed",
            "org.biojavax.bio.taxa.NCBITaxon",
            "NODERANK"
            );
    public static final ChangeType GENETICCODE = new ChangeType(
            "This taxon's genetic code has changed",
            "org.biojavax.bio.taxa.NCBITaxon",
            "GENETICCODE"
            );
    public static final ChangeType MITOGENETICCODE = new ChangeType(
            "This taxon's mito genetic code has changed",
            "org.biojavax.bio.taxa.NCBITaxon",
            "MITOGENETICCODE"
            );
    public static final ChangeType LEFTVALUE = new ChangeType(
            "This taxon's left value has changed",
            "org.biojavax.bio.taxa.NCBITaxon",
            "LEFTVALUE"
            );
    public static final ChangeType RIGHTVALUE = new ChangeType(
            "This taxon's right value has changed",
            "org.biojavax.bio.taxa.NCBITaxon",
            "RIGHTVALUE"
            );
    
	public static final ChangeType HIDDEN = new ChangeType(
			"This taxon's visibility in genbank hierarchy has changed",
			"org.biojavax.bio.taxa.NCBITaxon",
			"HIDDEN"
	);
	
    /**
     * Use this to define scientific names for things. There should
     * usually only be one scientific name for an organism.
     */
    public static final String SCIENTIFIC = "scientific name";
    
    /**
     * Use this to define common names for things. There can be as many
     * common names as you like.
     */
    public static final String COMMON = "common name";
    
    /**
     * Use this to define common names for things. There can be as many
     * synonyms as you like.
     */
    public static final String SYNONYM = "synonym";
    
    /**
     * Use this to define acronyms for things. There can be as many
     * acronyms as you like.
     */
    public static final String ACRONYM = "acronym";
    
    /**
     * Use this to define equivalent names for things. There can be as many
     * equivalent names as you like.
     */
    public static final String EQUIVALENT = "equivalent name";
    
    /**
     * Returns all the name classes available for a taxon. These are Strings.
     * @return a set of name classes, or the empty set if there are none.
     */
    public Set getNameClasses();
    
    /**
     * Returns all the names available for a taxon in a given class.
     * These are Strings.
     * @param nameClass the name class to retrieve names from.
     * @return a set of names, or the empty set if there are none.
     * @throws IllegalArgumentException if the name is null.
     */
    public Set getNames(String nameClass) throws IllegalArgumentException;
    
    /**
     * Adds the name to this taxon in the given name class. Neither can be null.
     * @param nameClass the class to add the name in.
     * @param name the name to add.
     * @throws ChangeVetoException in case of objections.
     * @throws IllegalArgumentException if the name is null.
     */
    public void addName(String nameClass, String name) throws IllegalArgumentException,ChangeVetoException;
    
    /**
     * Removes the name from the given name class. Neither can be null.
     * @return True if the name was found and removed, false otherwise.
     * @param nameClass the class to remove the name from.
     * @param name the name to remove.
     * @throws ChangeVetoException in case of objections.
     * @throws IllegalArgumentException if the name is null.
     */
    public boolean removeName(String nameClass, String name) throws IllegalArgumentException,ChangeVetoException;
    
    /**
     * Tests for the presence of a name in a given class. Neither can be null.
     * @param nameClass the class to look the name up in.
     * @param name the name to text for existence of.
     * @return True if the name exists in that class, false otherwise.
     * @throws IllegalArgumentException if the name is null.
     */
    public boolean containsName(String nameClass, String name) throws IllegalArgumentException;
    
    /**
     * Returns the parent NCBI taxon ID, if known.
     * @return Value of property parent.
     */
    public Integer getParentNCBITaxID();
    
    /**
     * Sets the parent NCBI taxon ID. May be null if not known.
     * @param parent New value of property parent.
     * @throws ChangeVetoException in case of objections.
     */
    public void setParentNCBITaxID(Integer parent) throws ChangeVetoException;
    
    /**
     * Gets the NCBI taxon ID. This is never null and is immutable, as otherwise
     * we would have no way of distinguishing between various taxa. It should
     * be set by the constructor of an implementation.
     * @return Value of property NCBITaxID.
     */
    public int getNCBITaxID();
    
    /**
     * Gets the node rank of this taxon. May be null.
     * @return Value of property nodeRank.
     */
    public String getNodeRank();
    
    /**
     * Sets the node rank of this taxon. May be null, in which case it is unset.
     * @param nodeRank New value of property nodeRank.
     * @throws ChangeVetoException in case of objections.
     */
    public void setNodeRank(String nodeRank) throws ChangeVetoException;
    
    /**
     * Returns the genetic code of this taxon, which may be null if not known.
     * @return Value of property geneticCode.
     */
    public Integer getGeneticCode();
    
    /**
     * Sets the genetic code of this taxon, which may be null, which will unset it.
     * @param geneticCode New value of property geneticCode.
     * @throws ChangeVetoException in case of objections.
     */
    public void setGeneticCode(Integer geneticCode) throws ChangeVetoException;
    
    /**
     * Returns the mitochondrial genetic code of this taxon, which may be null 
     * if not known.
     * @return Value of property mitoGeneticCode.
     */
    public Integer getMitoGeneticCode();
    
    /**
     * Sets the mitochondrial genetic code of this taxon, which may be null, 
     * which will unset it.
     * @param mitoGeneticCode New value of property mitoGeneticCode.
     * @throws ChangeVetoException in case of objections.
     */
    public void setMitoGeneticCode(Integer mitoGeneticCode) throws ChangeVetoException;
    
    /**
     * Gets the left value. May be null.
     * @return Value of property leftValue.
     */
    public Integer getLeftValue();
    
    /**
     * Sets the left value. May be null.
     * @param leftValue New value of property leftValue.
     * @throws ChangeVetoException in case of objections.
     */
    public void setLeftValue(Integer leftValue) throws ChangeVetoException;
    
    /**
     * Gets the right value. May be null.
     * @return Value of property rightValue.
     */
    public Integer getRightValue();
    
    /**
     * Sets the right value. May be null.
     * @param rightValue New value of property rightValue.
     * @throws ChangeVetoException in case of objections.
     */
    public void setRightValue(Integer rightValue) throws ChangeVetoException;

    /**
     * Returns the name of this taxon entry in the form:
     *   scientific (common)
     * or if there is no common name:
     *   scientific
     * @return the display name as described above.
     */
    public String getDisplayName();

    /**
     * Returns the taxonomy hierarchy of this taxon entry in the form:
     *   least specific; more specific; ...; most specific.
     * It follows the chain up the tree as far as it can, and will stop
     * at the first one it comes to that returns null for getParentNCBITaxID()
     * @return the display name as described above.
     */
    public String getNameHierarchy();
    
	/**
	 * used in getNameHierarchy() to determine whether this taxonomy level is displayed
	 */
	public boolean isTaxonHidden();
	
	/**
	 * determines whether this taxonomy level is displayed in etNameHierarchy()
	 * @param isTaxonHidden - if true it is displayed
	 * @throws ChangeVetoException
	 */
    public void setTaxonHidden(final boolean isTaxonHidden) throws ChangeVetoException;
}

