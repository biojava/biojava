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

package org.biojava.nbio.core.sequence.features;

import java.util.List;
import java.util.Map;

import org.biojava.nbio.core.sequence.location.template.AbstractLocation;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;

/**
 * Interface class to handle describing arbitrary features. A feature can be found at multiple locations in a sequence such as
 * the surface of a protein where different sequence positions make up that feature. Ligand binding pocket is another example.
 * The location in its current form knows the start and stop position in a sequence and thus should contain knowledge about the
 * actual sequence.
 *
 * A feature can contain features to handle cases where a domain is a feature and the secondary structures covered by that domain
 * and other requirements for grouping.
 *
 * @author Scooter Willis <willishf at gmail dot com>
 * @author Paolo Pavan
 */
public interface FeatureInterface<S extends AbstractSequence<C>, C extends Compound> {

	/**
	 * Get the short description that can be used to describe the feature
	 * @return
	 */

	public String getShortDescription();

	/**
	 * Set the short description that can be used to describe the feature
	 * @param shortDescription
	 */

	public void setShortDescription(String shortDescription);

	  /**
	 * Get the description that can be used to describe the feature
	 * @return
	 */

	public String getDescription();


	  /**
	 * Set the description that can be used to describe the feature
	 * @return
	 */

	public void setDescription(String description);

		/**
	 * The location(s) of this feature where the location should contain a reference to parent and sequence etc.
	 * <p>
	 * The location may be complicated, or simply a range.
	 * The annotation is assumed to apply to all the region contained
	 * within the location.
	 *
	 * @return a Location anchoring this feature
	 */
	public AbstractLocation getLocations();

		/**
	 * The new location for this feature.
	 * <p>
	 * The location may be complicated or simply a range. The annotation is
	 * assumed to apply to the entire region contained within the location.
	 * Any values returned from methods that rely on the old location must
	 * not be affected.
	 *
	 * @param loc the new Location for this feature
	 *
	 */
	public void setLocation(AbstractLocation loc);

		/**
	 * The type of the feature.
	 *
	 * @return the type of this sequence
	 */
	public String getType();

	/**
	 * Change the type of this feature.
	 *
	 * @param type  new type String
	 *
	 */
	public void setType(String type);


		/**
	 * The source of the feature. This may be a program or process.
	 *
	 * @return the source, or generator
	 */
	public String getSource();

	/**
	 * Change the source of the FeatureInterface.
	 *
	 * @param source the new source String
	 *
	 */
	public void setSource(String source);

	/**
	 * Set the parent feature
	 * @param feature
	 */

	public void setParentFeature(FeatureInterface<S, C> feature);

	/**
	 * Get the parent feature
	 * @return
	 */

	public FeatureInterface<S, C> getParentFeature();

	/**
	 * Get the features contained by this feature
	 * @return
	 */

	public List<FeatureInterface<S, C>> getChildrenFeatures();

	/**
	 * Set the children features
	 * @param features
	 */

	public void setChildrenFeatures(List<FeatureInterface<S, C>> features);


		/**
	 * @return the userObject
	 */
	public Object getUserObject();

	/**
	 * @param userObject the userObject to set
	 */
	public void setUserObject(Object userObject);

    /* new feature interface methods for accessing qualifers and the new map 
     * see abstract feature for implementation 
     */
    /**
     * map implementation to store qualifiers where only qualifier hold its key and value pair
     * @return
     */
    public GenBankQualifierMap getQualifierMap();
    /**
     * get all qualifiers of this feature
     * @return
     */
    public Qualifier[] getQualifiers();
	/**
	 * overwrite qualifiermap
	 * @param qualifierMap
	 */
	public void setQualifierMap(GenBankQualifierMap qualifierMap);
	/**
	 * overwrite qualifiers 
	 * @param qualifiers
	 */
	public void setQualifiers(Qualifier[] qualifiers);
	/**
	 * overwrite this qualifier 
	 * @param qualifiers
	 */
	public void setQualifier(Qualifier q);
	/**
	 * add this qualifier
	 * @param qualifier
	 */
	public void addQualifier(Qualifier qualifier);
	/**
	 * add a bunch of qualifiers  
	 * @param qa
	 */
	public void addQualifiers(Qualifier[] qa);
    public Qualifier getQualifierByName(String qName);
    public Qualifier getFirstQualifierByValue(String value);
    public Qualifier[] getQualifiersByValue(String value);
	//DBreferenceInfo 
	/**
	 * returns the dbreferenceinfo of this feature, which can contain lots of 
	 * entries 
	 * * @return
	 */
    public DBReferenceInfo getAllDatabaseReferenceInfos();
	/**
	 * returns all databases of this feature in a string[]
	 * @return
	 */
	public String[] getAllDatabases();
	/**
	 * returns all sequence database references for all databases in a string[]
	 * for this feature 
	 * @return
	 */
	public String[] getAllDatabaseReferences();
	/**
	 * get database reference info #i as new DBReferenceInfo
	 * @param i
	 * @return
	 */
	public DBReferenceInfo getDatabaseReferenceInfo(int i); 
	/**
	 * get database #i 
	 * @param i
	 * @return
	 */
	public String getDatabase(int i);
	/**
	 * get sequence database reference #i
	 * @param i
	 * @return
	 */
	public String getDatabaseReference(int i);
	/**
	 * convenience method to point out that there are several
	 * @return
	 */
	public DBReferenceInfo getFirstDatabaseReferenceInfo();
	public String getFirstDatabaseReference();
	public String getFirstDatabase();
	public String getDatabaseReference(String database, int i);
	public String[] getAllDatabaseReferences(String database);
	public String getFirstDatabaseReference(String database);
	public void setDatabaseReferenceInfo(DBReferenceInfo dbRefI);
	public void addDatabaseReferenceInfo(DBReferenceInfo dbRefI);
	//old stuff to be removed
	/**
	 * 
	 * @param str
	 * @param q
	 * Deprecated use addQualifier(Qualifier q)
	 */
	@Deprecated
	public void addQualifier(String str, Qualifier q);

}
