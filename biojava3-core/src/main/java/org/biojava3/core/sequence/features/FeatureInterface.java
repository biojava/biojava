/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.features;

import java.util.List;
import org.biojava3.core.sequence.location.SequenceLocation;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

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
    public SequenceLocation<S, C> getLocations();

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
    public void setLocation(SequenceLocation<S, C> loc);

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

}
