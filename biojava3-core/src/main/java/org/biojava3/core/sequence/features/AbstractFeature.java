/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.features;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import org.biojava3.core.sequence.location.SequenceLocation;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public abstract class AbstractFeature implements FeatureInterface {
    List<FeatureInterface> childrenFeatures = new ArrayList<FeatureInterface>();
    FeatureInterface parentFeature;
    SequenceLocation sequenceLocation;
    String type = "";
    String source = "";
    private String description = "";
    private String shortDescription = "";
    private Object userObject = null;

    public AbstractFeature(String type,String source){
        this.type = type;
        this.source = source;
    }

    @Override
    public SequenceLocation getLocations() {
        return sequenceLocation;
    }

    @Override
    public void setLocation(SequenceLocation loc) {
        sequenceLocation = loc;
    }

    @Override
    public String getType() {
        return type;
    }

    @Override
    public void setType(String type) {
        this.type = type;
    }

    @Override
    public String getSource() {
        return source;
    }

    @Override
    public void setSource(String source) {
        this.source = source;
    }

    @Override
    public void setParentFeature(FeatureInterface feature) {
        parentFeature = feature;
    }

    @Override
    public FeatureInterface getParentFeature() {
       return parentFeature;
    }

    @Override
    public List<FeatureInterface> getChildrenFeatures() {
        return childrenFeatures;
    }

    @Override
    public void setChildrenFeatures(List<FeatureInterface> features) {
        childrenFeatures = features;

    }

    /**
     * @return the description
     */
    public String getDescription() {
        return description;
    }

    /**
     * @param description the description to set
     */
    public void setDescription(String description) {
        this.description = description;
    }

    /**
     * @return the shortDescription
     */
    public String getShortDescription() {
        return shortDescription;
    }

    /**
     * @param shortDescription the shortDescription to set
     */
    public void setShortDescription(String shortDescription) {
        this.shortDescription = shortDescription;
    }

    /**
     * Sort features by start position and then longest length. When features are added
     * having them sorted by start position and then longest length helps on the layout
     * of overlapping features so they are delivered in a proper order.
     */

    static public final Comparator<FeatureInterface> LOCATION_LENGTH = new Comparator<FeatureInterface>() {

        public int compare(FeatureInterface e1, FeatureInterface e2) {
            double v1 = e1.getLocations().getStart().getPosition();
            double v2 = e2.getLocations().getStart().getPosition();
            if (v1 < v2) {
                return -1;
            } else if (v1 > v2) {
                return 1;
            } else {
                double end1 = e1.getLocations().getEnd().getPosition();
                double end2 = e2.getLocations().getEnd().getPosition();
                if(end1 > end2)
                    return -1;
                else if(end1 < end2)
                    return 1;
                else
                return 0;
            }

        }
    };

     /**
     * Sort features by length. //TODO need to handle cases where features have multiple locations, strand etc
     *
     */

    static public final Comparator<FeatureInterface> LENGTH = new Comparator<FeatureInterface>() {

        public int compare(FeatureInterface e1, FeatureInterface e2) {
            double v1 = Math.abs(e1.getLocations().getEnd().getPosition()- e1.getLocations().getStart().getPosition());
            double v2 = Math.abs(e2.getLocations().getEnd().getPosition() -  e2.getLocations().getStart().getPosition());
            if (v1 < v2) {
                return -1;
            } else if (v1 > v2) {
                return 1;
            } else {
                return 0;
            }

        }
    };

    /**
     * @return the userObject
     */
    public Object getUserObject() {
        return userObject;
    }

    /**
     * @param userObject the userObject to set
     */
    public void setUserObject(Object userObject) {
        this.userObject = userObject;
    }

}
