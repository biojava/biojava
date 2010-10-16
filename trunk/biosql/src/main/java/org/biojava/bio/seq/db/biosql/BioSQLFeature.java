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

package org.biojava.bio.seq.db.biosql;

import java.sql.SQLException;
import java.util.Iterator;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.RealizingFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.ontology.OntoTools;
import org.biojava.ontology.Term;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * @author Thomas Down
 * @author Simon Foote (modifications for schema version 1.0)
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 */
class BioSQLFeature implements Feature, RealizingFeatureHolder {
    private Annotation _annotation;
    private int id;

    // Feature stuff

    private String type;
    private String source;
    private Location location;

    // Relationship to sequences

    private int parentID = -1;
    private final BioSQLSequenceI sequence;

    // Children

    private SimpleFeatureHolder childFeatures;

    BioSQLFeature(Sequence seq,
		  Feature.Template templ)
	throws IllegalAlphabetException
    {
	this.type = templ.type;
	this.source = templ.source;
	this.location = templ.location;

	this.sequence = (BioSQLSequenceI) seq;

	_annotation = templ.annotation;
    }

    BioSQLFeature(Sequence seq,
		  FeatureHolder parent,
		  Feature.Template templ)
	throws IllegalAlphabetException
    {
	this(seq, templ);
	if (parent instanceof BioSQLFeature) {
	    parentID = ((BioSQLFeature) parent)._getInternalID();
	} else {
	    parentID = -1;
	}
    }

    public void hintChildFree() {
	if (childFeatures == null) {
	    childFeatures = new SimpleFeatureHolder();
	}
    }

    public void setParentID(int i) {
	this.parentID = i;
    }

    public void setType(String s)
        throws ChangeVetoException
    {
	BioSQLFeatureChangeHub featureHub = sequence.getSequenceDB().getFeatureChangeHub();
	ChangeEvent cev = new ChangeEvent(this, Feature.TYPE, getType(), s);
	synchronized (featureHub) {
	    featureHub.firePreChange(cev);
	    try {
		((BioSQLSequenceI) getSequence()).getSequenceDB().getFeaturesSQL().setFeatureType(id, s);
	    } catch (SQLException ex) {
		throw new BioRuntimeException("Error updating feature in database", ex);
	    }
	    this.type = s;
	    featureHub.firePostChange(cev);
	}
    }

    public String getType() {
	return type;
    }
    
    public Term getTypeTerm() {
        return OntoTools.ANY;
    }
    
    public void setTypeTerm(Term source) 
        throws ChangeVetoException
    {
        throw new ChangeVetoException("FIXME");
    }

    public void setSource(String s)
        throws ChangeVetoException
    {
        BioSQLFeatureChangeHub featureHub = sequence.getSequenceDB().getFeatureChangeHub();
	ChangeEvent cev = new ChangeEvent(this, Feature.SOURCE, getSource(), s);
	synchronized (featureHub) {
            featureHub.firePreChange(cev);
	    try {
		((BioSQLSequenceI) getSequence()).getSequenceDB().getFeaturesSQL().setFeatureSource(id, s);
	    } catch (SQLException ex) {
		throw new BioRuntimeException("Error updating feature in database", ex);
	    }
	    this.source = s;
            featureHub.firePostChange(cev);
	}
    }

    public String getSource() {
	return source;
    }
    
    public Term getSourceTerm() {
        return OntoTools.ANY;
    }
    
    public void setSourceTerm(Term source) 
        throws ChangeVetoException
    {
        throw new ChangeVetoException("FIXME");
    }

    public void setLocation(Location l)
        throws ChangeVetoException
    {
        BioSQLFeatureChangeHub featureHub = sequence.getSequenceDB().getFeatureChangeHub();
	ChangeEvent cev = new ChangeEvent(this, Feature.LOCATION, getLocation(), l);
	synchronized (featureHub) {
	    featureHub.firePreChange(cev);
	    try {
		((BioSQLSequenceI) getSequence()).getSequenceDB().getFeaturesSQL().setFeatureLocation(id, l, StrandedFeature.UNKNOWN);
	    } catch (SQLException ex) {
		throw new BioRuntimeException("Error updating feature in database", ex);
	    }
	    this.location = l;
            featureHub.firePostChange(cev);
	}
    }

    public Location getLocation() {
	return location;
    }

    public FeatureHolder getParent() {
	if (parentID == -1) {
	    return sequence;
	} else {
	    return sequence.getSequenceDB().getFeatureByID(parentID);
	}
    }

    public Sequence getSequence() {
	return sequence;
    }

    public void _setAnnotation(Annotation a) {
	_annotation = a;
    }

    public Feature realizeFeature(FeatureHolder fh, Feature.Template templ)
        throws BioException
    {
	try {
	    RealizingFeatureHolder rfh = (RealizingFeatureHolder) getParent();
	    return rfh.realizeFeature(fh, templ);
	} catch (ClassCastException ex) {
	    throw new BioException("Couldn't propagate feature creation request.");
	}
    }

    public Feature createFeature(Feature.Template templ)
        throws BioException, ChangeVetoException
    {
	Feature f = realizeFeature(this, templ);
	
        BioSQLFeatureChangeHub featureHub = sequence.getSequenceDB().getFeatureChangeHub();
	ChangeEvent cev = new ChangeEvent(this, FeatureHolder.FEATURES, f, null);
	synchronized (featureHub) {
            featureHub.firePreChange(cev);
	    getFeatures().addFeature(f);
	    ((BioSQLSequenceI) getSequence()).persistFeature(f, id);
            featureHub.firePostChange(cev);
	}
	return f;
    }

    public void removeFeature(Feature f)
        throws ChangeVetoException
    {
        BioSQLFeatureChangeHub featureHub = sequence.getSequenceDB().getFeatureChangeHub();
	ChangeEvent cev = new ChangeEvent(this, FeatureHolder.FEATURES, null, f);
	synchronized (featureHub) {
            featureHub.firePreChange(cev);
	    getFeatures().removeFeature(f);
	    ((BioSQLSequenceI) getSequence()).getSequenceDB().getFeaturesSQL().removeFeature((BioSQLFeature) f);
	    featureHub.firePostChange(cev);
	}
    }

    public Annotation getAnnotation() {
	return _annotation;
    }

    public void _setInternalID(int i) {
	this.id = i;
    }

    public int _getInternalID() {
	return id;
    }

    public synchronized void _addFeature(Feature f) 
        throws ChangeVetoException
    {
	if (childFeatures == null) {
	    childFeatures = new SimpleFeatureHolder();
	}
	childFeatures.addFeature(f);
    }

    protected void fillTemplate(Feature.Template template) {
	template.source = source;
	template.type = type;
	template.location = location;
	template.annotation = _annotation;
    }

    public Feature.Template makeTemplate() {
	Feature.Template template = new Feature.Template();
	fillTemplate(template);
	return template;
    }

    public Iterator features() {
	return getFeatures().features();
    }

    public int countFeatures() {
	return getFeatures().countFeatures();
    }

    public boolean containsFeature(Feature f) {
	return getFeatures().containsFeature(f);
    }

    public FeatureFilter getSchema() {
        return new FeatureFilter.ByParent(new FeatureFilter.ByFeature(this));
    }
    
    public FeatureHolder filter(FeatureFilter ff) {
        FeatureFilter childFilter = new FeatureFilter.And(new FeatureFilter.ContainedByLocation(getLocation()),
                                                          new FeatureFilter.Not(FeatureFilter.top_level));
                                                          
        if (FilterUtils.areDisjoint(ff, childFilter)) {
            return FeatureHolder.EMPTY_FEATURE_HOLDER;
        } else {
            return getFeatures().filter(ff);
        }
    }
    
    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
        FeatureFilter childFilter = new FeatureFilter.ContainedByLocation(getLocation());
        if (FilterUtils.areDisjoint(ff, childFilter)) {
            return FeatureHolder.EMPTY_FEATURE_HOLDER;
        } else {
            return getFeatures().filter(ff, recurse);
        }
    }

    private class FeatureReceiver extends BioSQLFeatureReceiver {
	FeatureReceiver() {
	    super(sequence);
	}

	protected void deliverTopLevelFeature(Feature f)
	    throws ParseException, ChangeVetoException
	{
	    childFeatures.addFeature(f);
	}
    }

    private synchronized SimpleFeatureHolder getFeatures() {
	if (childFeatures == null) {
	    try {
		BioSQLSequenceI seqi = (BioSQLSequenceI) sequence;
		childFeatures = new SimpleFeatureHolder();
		FeaturesSQL adaptor = seqi.getSequenceDB().getFeaturesSQL();
		adaptor.retrieveFeatures(seqi.getBioEntryID(),
					 new FeatureReceiver(),
					 null,
					 id,
					 -1);
	    } catch (SQLException ex) {
		throw new BioRuntimeException("SQL error while reading features", ex);
	    } catch (BioException ex) {
		throw new BioRuntimeException(ex);
	    } 
	}

	return childFeatures;
    }

    public SymbolList getSymbols() {
	return getLocation().symbols(getSequence());
    }

    public int hashCode() {
	return makeTemplate().hashCode();
    }

//      public boolean equals(Object o) {
//  	if (! (o instanceof Feature)) {
//  	    return false;
//  	}

//  	Feature fo = (Feature) o;
//  	if (! fo.getSequence().equals(getSequence())) 
//  	    return false;
    
//  	return makeTemplate().equals(fo.makeTemplate());
//      }
    
    public void addChangeListener(ChangeListener cl) {
	addChangeListener(cl, ChangeType.UNKNOWN);
    }
    
    public void addChangeListener(ChangeListener cl, ChangeType ct) {
	sequence.getSequenceDB().getFeatureChangeHub().addListener(new Integer(id), cl, ct);
    }

    public void removeChangeListener(ChangeListener cl) {
	removeChangeListener(cl, ChangeType.UNKNOWN);
    }

    public void removeChangeListener(ChangeListener cl, ChangeType ct) {
	sequence.getSequenceDB().getFeatureChangeHub().removeListener(new Integer(id), cl, ct);
    }

    public boolean isUnchanging(ChangeType ct) {
	return false;
    }
} 
