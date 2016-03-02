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

package org.biojava.nbio.core.sequence.features;

import org.biojava.nbio.core.sequence.location.template.AbstractLocation;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * It is {@link DBReferenceInfo} which implements {@link FeatureInterface}. It allows to keep a dbReferenceInfo as a feature.
 * 
 * @author Jacek Grzebyta
 * @author Paolo Pavan
 * @param <S>
 * @param <C>
 */
public class FatureDbReferenceInfo<S extends AbstractSequence<C>, C extends Compound> extends DBReferenceInfo implements FeatureInterface<S,C> {
    
    private AbstractLocation location;
    private FeatureInterface<S,C> parentFeature;
    private List<FeatureInterface<S, C>> childrenFeatures = new ArrayList<FeatureInterface<S, C>>();
    private String description = "";
    private String shortDescription = "";
    private Object userObject;
    private GenBankQualifierMap qualifierMap = new GenBankQualifierMap();
    
    
    public FatureDbReferenceInfo(String database, String id) {
        super(database+":"+id);
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
    public AbstractLocation getLocations() {
        return location;
    }

    @Override
    public void setLocation(AbstractLocation loc) {
        location = loc;
    }

/*
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
*/
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

    public void setQualifierMap(GenBankQualifierMap qMap) {
        this.qualifierMap = qMap;
    }


    public void addQualifier(Qualifier qualifier) {
    	qualifierMap.add(qualifier);
    }

    public void addQualifiers(Qualifier[] qa) {
    	qualifierMap.addQualifiers(qa);
    }
    
    public void setQualifiers(Qualifier[] qa) {
    	qualifierMap.set(qa);
    }
    
	public Qualifier[] getQualifiers() {
		// TODO Auto-generated method stub
		return null;
	}

	
	public GenBankQualifierMap getQualifierMap() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getType() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setType(String type) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String getSource() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setSource(String source) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Qualifier getQualifierByName(String qName) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Qualifier getFirstQualifierByValue(String value) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Qualifier[] getQualifiersByValue(String value) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void addQualifier(String str, Qualifier q) {
		// TODO Auto-generated method stub
		
	}
}