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

package org.biojava.bio.seq.distributed;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.AbstractSequenceDB;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Sequence database from the meta-DAS system.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.2
 *
 * Once you've made one of these and populated it with a few DistDataSource instances,
 * you should be able to prety much forget about it and use it directly as a normal
 * SequenceDB implementation.
 *
 * DataSources can be added and removed while the object is live. 
 */

public class DistributedSequenceDB extends AbstractSequenceDB implements SequenceDB {
    public static final ChangeType DATASOURCE = new ChangeType(
	    "Data sources have changes in a Distributed Sequence DB",
	    "org.biojava.bio.seq.distributed.DistributedSequenceDB",
	    "DATASOURCE",
	    ChangeType.UNKNOWN
    );

    public static final ChangeType DATASOURCE_SELECTION = new ChangeType(
	    "The set of available data sources has changes in a Distributed Sequence DB",
	    "org.biojava.bio.seq.distributed.DistributedSequenceDB",
	    "DATASOURCE_SELECTION",
	    DistributedSequenceDB.DATASOURCE
    );

    private Set datasources;
    private transient ChangeSupport changeSupport;

    protected boolean hasChangeSupport() {
	return (changeSupport != null);
    }

    protected ChangeSupport getChangeSupport() {
	if (changeSupport == null) {
	    changeSupport = new ChangeSupport();
	}
	return changeSupport;
    }

    {
	datasources = new HashSet();
    }

    /**
     * Get the current set of all currently registered data sources.
     *
     * <p>
     * The returned Set is totaly independant of any future changes made to the
     * distributed sequence DB.
     * </p>
     *
     * @return a new Set containing all DistDataSource instances registered at
     *    the time
     */ 
    public Set getDataSources() {
	return new HashSet(datasources);
    }

    /**
     * Add a distributed data source.
     *
     * @param dds  the DistDataSource to add
     * @throws ChangeVetoException if the data source could not be added
     */
    public void addDataSource(DistDataSource dds) 
        throws ChangeVetoException
    {
	if (datasources.contains(dds)) {
	    return;
	}

	if (hasChangeSupport()) {
	    ChangeSupport cs = getChangeSupport();
	    synchronized (cs) {
		ChangeEvent cev = new ChangeEvent(this,
						  DATASOURCE_SELECTION,
						  dds,
						  null);
		cs.firePreChangeEvent(cev);
		_addDataSource(dds);
		cs.firePostChangeEvent(cev);
	    }
	} else {
	    _addDataSource(dds);
	}
    }

    private void _addDataSource(DistDataSource dds) {
	datasources.add(dds);
    }

    /**
     * Remove a distributed data source.
     *
     * @param dds  the DistributedDataSource to remove
     * @throws ChangeVetoException  if the data source could not be removed
     */
    public void removeDataSource(DistDataSource dds) 
        throws ChangeVetoException
    {
	if (!datasources.contains(dds)) {
	    throw new ChangeVetoException("That datasource isn't currently installed");
	}

	if (hasChangeSupport()) {
	    ChangeSupport cs = getChangeSupport();
	    synchronized (cs) {
		ChangeEvent cev = new ChangeEvent(this,
						  DATASOURCE_SELECTION,
						  null,
						  dds);
		cs.firePreChangeEvent(cev);
		_removeDataSource(dds);
		cs.firePostChangeEvent(cev);
	    }
	} else {
	    _removeDataSource(dds);
	}
    }

    private void _removeDataSource(DistDataSource dds) {
	datasources.remove(dds);
    }

    public String getName() {
	return "<unknown meta-das>";
    }

    public void addSequence(Sequence seq)
        throws ChangeVetoException
    {
	throw new ChangeVetoException("Can't add sequences to meta-das");
    }

    public void removeSequence(String id)
        throws ChangeVetoException
    {
	throw new ChangeVetoException("Can't add sequences to meta-das");
    }

    public Sequence getSequence(String id)
        throws IllegalIDException, BioException
    {
	Set featureSources = new HashSet();
	DistDataSource seqSource = null;

	for (Iterator i = datasources.iterator(); i.hasNext(); ) {
	    DistDataSource dds = (DistDataSource) i.next();
	    if (dds.hasSequence(id) && seqSource == null) {
		seqSource = dds;
	    }

	    if (dds.hasFeatures(id)) {
		featureSources.add(dds);
	    }
	}

	if (seqSource == null) {
	    throw new IllegalIDException("No sequence source for ID: " + id);
	}

	return new DistributedSequence(id, this, seqSource, featureSources);
    }

    public Set ids() {
	Set ids = new HashSet();
	for (Iterator i = datasources.iterator(); i.hasNext(); ) {
	    DistDataSource dds = (DistDataSource) i.next();
	    try {
		ids.addAll(dds.ids(true));
	    } catch (BioException ex) {
	    }
	}
	return ids;
    }

}
