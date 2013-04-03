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
 * Created on 01-21-2010
 */

package org.biojava3.core.sequence.features;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import org.biojava3.core.sequence.location.SequenceLocation;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 * A feature is currently any descriptive item that can be associated with a sequence position(s)
 * A feature has a type and a source which is currently a string to allow flexibility for the user
 * Ideally well defined features should have a class to describe attributes of that feature
 * @author Scooter Willis <willishf at gmail dot com>
 */
public abstract class AbstractFeature<S extends AbstractSequence<C>, C extends Compound>
        implements FeatureInterface<S, C> {
    List<FeatureInterface<S, C>> childrenFeatures = new ArrayList<FeatureInterface<S, C>>();
    FeatureInterface<S, C> parentFeature;
    SequenceLocation<S, C> sequenceLocation;
    String type = "";
    String source = "";
    private String description = "";
    private String shortDescription = "";
    private Object userObject = null;

    /**
     * A feature has a type and a source
     * @param type
     * @param source
     */
    public AbstractFeature(String type,String source){
        this.type = type;
        this.source = source;
    }

    /**
     * A feature could be a single sequence position like a mutation or a post translational modification of an amino acid.
     * It could also be the docking interface of N number of amino acids on the surface. The location wold then be a collection
     * of sequence positions instead of a single sequence position or the begin and end of a sequence seqment.
     * @return
     */

    @Override
    public SequenceLocation<S, C> getLocations() {
        return sequenceLocation;
    }

    /**
     *  A feature could be a single sequence position like a mutation or a post translational modification of an amino acid.
     * It could also be the docking interface of N number of amino acids on the surface. The location wold then be a collection
     * of sequence positions instead of a single sequence position or the begin and end of a sequence seqment.
     * @param loc
     */
    @Override
    public void setLocation(SequenceLocation<S, C> loc) {
        sequenceLocation = loc;
    }

    /**
     * The feature type
     * @return
     */
    @Override
    public String getType() {
        return type;
    }

    /**
     * Set the feature type
     * @param type
     */
    @Override
    public void setType(String type) {
        this.type = type;
    }

    /**
     * The feature source
     * @return
     */

    @Override
    public String getSource() {
        return source;
    }

    /**
     * Set the feature source
     * @param source
     */
    @Override
    public void setSource(String source) {
        this.source = source;
    }

    /**
     * A feature can be the child or contained by a parent feature. An example is a Helix feature could contain
     * children features. A PFAM domain could contain secondary structures.
     * @param feature
     */
    @Override
    public void setParentFeature(FeatureInterface<S, C> feature) {
        parentFeature = feature;
    }

    /**
     * Get the parent Feature
     * @return
     */
    @Override
    public FeatureInterface<S, C> getParentFeature() {
       return parentFeature;
    }

    /**
     * Get the children features
     * @return
     */
    @Override
    public List<FeatureInterface<S, C>> getChildrenFeatures() {
        return childrenFeatures;
    }

    /**
     * Set the children features
     * @param features
     */
    @Override
    public void setChildrenFeatures(List<FeatureInterface<S, C>> features) {
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

    public static final Comparator<FeatureInterface<?, ?>> LOCATION_LENGTH = new Comparator<FeatureInterface<?, ?>>() {

        public int compare(FeatureInterface<?, ?> e1, FeatureInterface<?, ?> e2) {
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

    static public final Comparator<FeatureInterface<?, ?>> LENGTH = new Comparator<FeatureInterface<?, ?>>() {

        public int compare(FeatureInterface<?, ?> e1, FeatureInterface<?, ?> e2) {
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
     * Allow the user to associate an object with the feature. This way if a feature which is displayed in a GUI
     * is clicked on the application can then get a user defined object associated with the feature.
     * @param userObject the userObject to set
     */
    public void setUserObject(Object userObject) {
        this.userObject = userObject;
    }

}
