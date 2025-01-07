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
 * <p>
 * A feature can contain features to handle cases where a domain is a feature and the secondary structures covered by that domain
 * and other requirements for grouping.
 *
 * @author Scooter Willis 
 * @author Paolo Pavan
 */
public interface FeatureInterface<S extends AbstractSequence<C>, C extends Compound> {

	/**
	 * Get the short description that can be used to describe the feature
	 * @return
	 */
	String getShortDescription();

	/**
	 * Set the short description that can be used to describe the feature
	 * @param shortDescription
	 */
	void setShortDescription(String shortDescription);

	 /**
	 * Get the description that can be used to describe the feature
	 * @return
	 */
	String getDescription();

	 /**
	 * Set the description that can be used to describe the feature
	 */
	void setDescription(String description);

	/**
	 * The location(s) of this feature where the location should contain a reference to parent and sequence etc.
	 * <p>
	 * The location may be complicated, or simply a range.
	 * The annotation is assumed to apply to all the region contained
	 * within the location.
	 *
	 * @return a Location anchoring this feature
	 */
	AbstractLocation getLocations();

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
	void setLocation(AbstractLocation loc);

	/**
	 * The type of the feature.
	 *
	 * @return the type of this sequence
	 */
	String getType();

	/**
	 * Change the type of this feature.
	 *
	 * @param type  new type String
	 *
	 */
	void setType(String type);


	/**
	 * The source of the feature. This may be a program or process.
	 *
	 * @return the source, or generator
	 */
	String getSource();

	/**
	 * Change the source of the FeatureInterface.
	 *
	 * @param source the new source String
	 *
	 */
	void setSource(String source);

	/**
	 * Set the parent feature
	 * @param feature
	 */
	void setParentFeature(FeatureInterface<S, C> feature);

	/**
	 * Get the parent feature
	 * @return
	 */
	FeatureInterface<S, C> getParentFeature();

	/**
	 * Get the features contained by this feature
	 * @return
	 */
	List<FeatureInterface<S, C>> getChildrenFeatures();

	/**
	 * Set the children features
	 * @param features
	 */
	void setChildrenFeatures(List<FeatureInterface<S, C>> features);


	/**
	 * @return the userObject
	 */
	Object getUserObject();

	/**
	 * @param userObject the userObject to set
	 */
	void setUserObject(Object userObject);

	/**
	 * Get the qualifiers for this feature
	 * @return
	 */

	Map<String, List<Qualifier>> getQualifiers();

	/**
	 * Set the qualifiers
	 * @param qualifiers
	 */
	void setQualifiers(Map<String, List<Qualifier>> qualifiers);

	/**
	 * Add a qualifier
	 * @param qualifier
	 */
	void addQualifier(String key, Qualifier qualifier);

}
