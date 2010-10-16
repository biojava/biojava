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

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.MergeFeatureHolder;
import org.biojava.bio.seq.RealizingFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.cache.CacheReference;

/**
 * Top-level SeqFeature set for a BioEntry
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */

class BioSQLTiledFeatures implements FeatureHolder, RealizingFeatureHolder
{
    private Sequence seq;
    private BioSQLSequenceDB seqDB;
    private int bioentry_id;

    private Location[]            tileLocations;
    private FeatureTile[]         tileFeatures;
    private SimpleFeatureHolder   overlappingFeatures;
    private MergeFeatureHolder    allFeatures;
	
    BioSQLTiledFeatures(Sequence seq,
			BioSQLSequenceDB seqDB,
			int bioentry_id,
			int tileSize)
    {
	this.seq = seq;
	this.seqDB = seqDB;
	this.bioentry_id = bioentry_id;

	int numTiles = (int) Math.ceil((1.0 * seq.length()) / tileSize);
	tileLocations = new Location[numTiles];
	tileFeatures = new FeatureTile[numTiles];

	try {
	    allFeatures = new MergeFeatureHolder();
	    for (int t = 0; t < numTiles; ++t) {
		tileLocations[t] = new RangeLocation(1 + (t * tileSize),
						     Math.min((t + 1) * tileSize, seq.length()));
		tileFeatures[t] = new FeatureTile(t);
		allFeatures.addFeatureHolder(tileFeatures[t]);
	    }
	    
	    overlappingFeatures = new SimpleFeatureHolder();
	    allFeatures.addFeatureHolder(overlappingFeatures);
	} catch (ChangeVetoException ex) {
	    throw new BioError(ex);
	}
    }

    
    public FeatureFilter getSchema() {
        return FeatureFilter.top_level;
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

    public FeatureHolder filter(FeatureFilter ff) {
        return getFeatures().filter(ff);
    }
    
    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
        return getFeatures().filter(ff, recurse);
    }

    private void _addFeature(Feature f) 
        throws ChangeVetoException
    {
	for (int t = 0; t < tileLocations.length; ++t) {
	    if (tileLocations[t].contains(f.getLocation())) {
		tileFeatures[t].addFeature(f);
		return;
	    }
	}

	overlappingFeatures.addFeature(f);
    }

    public Feature createFeature(Feature.Template ft)
        throws ChangeVetoException, BioException
    {
	Feature f = realizeFeature(seq, ft);
	
	BioSQLEntryChangeHub entryHub = ((BioSQLSequenceI) seq).getSequenceDB().getEntryChangeHub();
	ChangeEvent cev = new ChangeEvent(seq, FeatureHolder.FEATURES, f);
	synchronized (entryHub) {
	    entryHub.firePreChange(cev);
	    seqDB.getFeaturesSQL().persistFeature(f, -1, bioentry_id); // No parent
	    _addFeature(f);
	    entryHub.firePostChange(cev);
	}

	return f;
    }

    public void removeFeature(Feature f)
        throws ChangeVetoException, BioException
    {
	FeatureHolder fh = overlappingFeatures;
	for (int t = 0; t < tileLocations.length; ++t) {
	    if (tileLocations[t].contains(f.getLocation())) {
		fh = tileFeatures[t];
	    }
	}

        if (!fh.containsFeature(f)) {
            throw new ChangeVetoException("Feature doesn't come from this sequence");
        }
        if (!(f instanceof BioSQLFeature)) {
            throw new ChangeVetoException("This isn't a normal BioSQL feature");
        }
        
        BioSQLEntryChangeHub entryHub = ((BioSQLSequenceI) seq).getSequenceDB().getEntryChangeHub();
	ChangeEvent cev = new ChangeEvent(seq, FeatureHolder.FEATURES, f);
	synchronized (entryHub) {
	    entryHub.firePreChange(cev);
	    seqDB.getFeaturesSQL().removeFeature((BioSQLFeature) f);
	    fh.removeFeature(f);
	    entryHub.firePostChange(cev);
        }
    }

    protected FeatureHolder getFeatures() {
	return allFeatures;
    }

    //
    // implements RealizingFeatureHolder
    //

    private BioSQLFeature _realizeFeature(FeatureHolder parent, Feature.Template templ)
        throws BioException
    {
	if (parent != seq && !seqDB.isHierarchySupported()) {
	    throw new BioException("This database doesn't support feature hierarchy.  Please create a seqfeature_relationship table");
	}

	if (templ instanceof StrandedFeature.Template && seq.getAlphabet() == DNATools.getDNA()) {
	    return new BioSQLStrandedFeature(seq, parent, (StrandedFeature.Template) templ);
	} else {
	    return new BioSQLFeature(seq, parent, templ);
	}
    }

    public Feature realizeFeature(FeatureHolder parent, Feature.Template templ)
        throws BioException
    {
	return _realizeFeature(parent, templ);
    }

    private class TileFeatureReceiver extends BioSQLFeatureReceiver {
	private final Location tileLocation;
	private final SimpleFeatureHolder tileFeatures;

	private TileFeatureReceiver(SimpleFeatureHolder tileFeatures,
				    Location tileLocation)
	{
	    super(seq);
	    this.tileLocation = tileLocation;
	    this.tileFeatures = tileFeatures;
	}

	protected void deliverTopLevelFeature(Feature f)
	    throws ParseException, ChangeVetoException
	{
	    if (LocationTools.contains(tileLocation, f.getLocation())) {
		tileFeatures.addFeature(f);
	    } else {
		if (!overlappingFeatures.containsFeature(f)) {
		    // System.err.println("Adding feature at " + f.getLocation() + " to overlaps");
		    overlappingFeatures.addFeature(f);
		}
	    }
	}
    }

    private class FeatureTile implements FeatureHolder {
	private int tileNumber;
	private CacheReference featuresRef;

	FeatureTile(int tileNumber) {
	    this.tileNumber = tileNumber;
	}

	public Iterator features() {
	    return getTileFeatures().features();
	}
	
	public int countFeatures() {
	    return getTileFeatures().countFeatures();
	}

	public boolean containsFeature(Feature f) {
	    return getTileFeatures().containsFeature(f);
	}
    
    public FeatureHolder filter(FeatureFilter ff) {
	    return getTileFeatures().filter(ff);
	}

	public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
	    return getTileFeatures().filter(ff, recurse);
	}

	private synchronized SimpleFeatureHolder getTileFeatures() {
	    if (featuresRef != null) {
		SimpleFeatureHolder fh = (SimpleFeatureHolder) featuresRef.get();
		if (fh == null) {
		    // System.err.println("*** Tile cache was cleared: " + tileNumber);
		} else {
		    return fh;
		}
	    }

	    try {
		// System.err.println("*** Fetching" + tileNumber);
		SimpleFeatureHolder features = new SimpleFeatureHolder();
		FeaturesSQL adaptor = seqDB.getFeaturesSQL();
		adaptor.retrieveFeatures(bioentry_id, 
					 new TileFeatureReceiver(features, tileLocations[tileNumber]),
					 tileLocations[tileNumber],
					 -1,
					 -1);
		featuresRef = seqDB.getTileCache().makeReference(features);
		return features;
	    } catch (SQLException ex) {
		throw new BioRuntimeException("SQL error while reading features", ex);
	    } catch (BioException ex) {
		throw new BioRuntimeException(ex);
	    } 
	}
    
    public FeatureFilter getSchema() {
        FeatureFilter tileFilter = new FeatureFilter.ContainedByLocation(tileLocations[tileNumber]);
        return new FeatureFilter.And(
                tileFilter,
                new FeatureFilter.OnlyDescendants(tileFilter)
        );
    }

	public void addFeature(Feature f) 
	    throws ChangeVetoException
	{
	    getTileFeatures().addFeature(f);
	}

	public void removeFeature(Feature f) 
	    throws ChangeVetoException
	{
	    getTileFeatures().removeFeature(f);
	}
	    
	public Feature createFeature(Feature.Template ft)
	    throws ChangeVetoException, BioException
	{
	    throw new ChangeVetoException();
	}

	// Not changeable
    
	public void addChangeListener(ChangeListener cl) {}
	public void addChangeListener(ChangeListener cl, ChangeType ct) {}
	public void removeChangeListener(ChangeListener cl) {}
	public void removeChangeListener(ChangeListener cl, ChangeType ct) {}
	public boolean isUnchanging(ChangeType ct) {
	    return true;
	}
    }

    
    public void addChangeListener(ChangeListener cl) {}
    public void addChangeListener(ChangeListener cl, ChangeType ct) {}
    public void removeChangeListener(ChangeListener cl) {}
    public void removeChangeListener(ChangeListener cl, ChangeType ct) {}
    public boolean isUnchanging(ChangeType ct) {
	return true;
    }
}
