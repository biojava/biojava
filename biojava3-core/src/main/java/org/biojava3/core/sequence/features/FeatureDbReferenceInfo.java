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
 * @author Jacek Grzebyta <github:jgrzebyta>
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 25-07-2014
 *
 */

package org.biojava3.core.sequence.features;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.biojava3.core.sequence.location.SequenceLocation;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 * It is {@link DBReferenceInfo} which implements {@link FeatureInterface}. It allows to keep a dbReferenceInfo as a feature.
 * 
 * @author Jacek Grzebyta
 * @param <S>
 * @param <C>
 */
public class FeatureDbReferenceInfo<S extends AbstractSequence<C>, C extends Compound> extends DBReferenceInfo implements FeatureInterface<S,C> {
    
    private SequenceLocation<S,C> location;
    private FeatureInterface<S,C> parentFeature;
    private List<FeatureInterface<S, C>> childrenFeatures = new ArrayList<FeatureInterface<S, C>>();
    private String description = "";
    private String shortDescription = "";
    private Object userObject;
    private HashMap<String,Qualifier> qualifiers = new HashMap<String,Qualifier>();
    
    
    public FeatureDbReferenceInfo(String database, String id) {
        super(database, id);
    }
    
    @Override
    public String getShortDescription() {
        return shortDescription;
    }

    @Override
    public void setShortDescription(String shortDescription) {
        this.shortDescription = shortDescription;
    }

    @Override
    public String getDescription() {
        return description;
    }

    @Override
    public void setDescription(String description) {
        this.description = description;
    }

	@Override
    public SequenceLocation<S, C> getLocations() {
        return location;
    }

    @Override
    public void setLocation(SequenceLocation<S, C> loc) {
        location = loc;
    }

    @Override
    public String getType() {
        return super.getDatabase();
    }

    @Override
    public void setType(String type) {
       super.setDatabase(type);
    }

    @Override
    public String getSource() {
        return super.getId();
    }

    @Override
    public void setSource(String source) {
        super.setId(source);
    }

    @Override
    public void setParentFeature(FeatureInterface<S,C> feature) {
        this.parentFeature = feature;
    }

    @Override
    public FeatureInterface<S,C> getParentFeature() {
        return this.parentFeature;
    }

    @Override
    public List<FeatureInterface<S, C>> getChildrenFeatures() {
        return this.childrenFeatures;
    }

    @Override
    public void setChildrenFeatures(List<FeatureInterface<S, C>> features) {
        this.childrenFeatures = features;
    }

    @Override
    public Object getUserObject() {
        return this.userObject;
    }

    @Override
    public void setUserObject(Object userObject) {
        this.userObject = userObject;
    }

    @Override
    public HashMap<String, Qualifier> getQualifiers() {
        return qualifiers;
    }

    @Override
    public void setQualifiers(HashMap<String, Qualifier> qualifiers) {
        this.qualifiers = qualifiers;
    }

    @Override
    public void addQualifier(String key, Qualifier qualifier) {
        if (qualifiers == null) {
            qualifiers = new HashMap<String, Qualifier>();
        }
        
        qualifiers.put(key, qualifier);
    }
    
}
