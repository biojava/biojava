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

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.SmallAnnotation;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Annotation keyed off a BioSQL seqfeature_qualifier_value table
 *
 * @author Thomas Down
 * @author Len Trigg
 * @author Richard Holland
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */

class BioSQLFeatureAnnotation
  implements Annotation
{
    private BioSQLSequenceDB seqDB;
    private int feature_id;
    private Annotation underlyingAnnotation;

    BioSQLFeatureAnnotation(BioSQLSequenceDB seqDB,
			    int feature_id)
    {
	this.seqDB = seqDB;
	this.feature_id = feature_id;
    }

    int getFeatureID() {
	return feature_id;
    }

    private void initAnnotations() {
      Connection conn = null;
	try {
	    conn = seqDB.getDataSource().getConnection();

	    PreparedStatement get_annotations = conn.prepareStatement("select term.name, seqfeature_qualifier_value.value " +
								      "  from term, seqfeature_qualifier_value " +
								      " where seqfeature_qualifier_value.seqfeature_id = ? and " +
								      "       term.term_id = seqfeature_qualifier_value.term_id");
	    get_annotations.setInt(1, feature_id);
	    ResultSet rs = get_annotations.executeQuery();
	    
            underlyingAnnotation = new SmallAnnotation(); // Initialize here in case there are no resultset entries.
	    while (rs.next()) {
		String key = rs.getString(1).trim();   // HACK due to stupid schema change
		String value = rs.getString(2);
		try {
		    initProperty(key, value);
		} catch (ChangeVetoException ex) {
		    rs.close();
                get_annotations.close();
                conn.close();
		    throw new BioError(ex);
		}
	    }
            rs.close();
	    get_annotations.close();
	    conn.close();
	} catch (SQLException ex) {
          if (conn!=null) try {conn.close();} catch (SQLException ex3) {}
	    throw new BioRuntimeException("Error fetching annotations", ex);
	}
    }

    public Object getProperty(Object key) {
	if (underlyingAnnotation == null) {
	    initAnnotations();
	}

	return underlyingAnnotation.getProperty(key);
    }

    void initProperty(Object key, Object value)
        throws ChangeVetoException
    {
        if (underlyingAnnotation == null) {
            underlyingAnnotation = new SmallAnnotation();
        }
        if (underlyingAnnotation.containsProperty(key)) {
            Object oldVal = underlyingAnnotation.getProperty(key);
            if (oldVal instanceof Collection) {
                ((Collection) oldVal).add(value);
            } else {
                List newList = new ArrayList();
                newList.add(oldVal);
                newList.add(value);
                underlyingAnnotation.setProperty(key, newList);
            }
        } else {
            underlyingAnnotation.setProperty(key, value);
        }
    }
    
    public void setProperty(Object key, Object value)
        throws ChangeVetoException
    {
        BioSQLFeatureAnnotationChangeHub featureAnnotationHub = seqDB.getFeatureAnnotationChangeHub();
        synchronized (featureAnnotationHub) {
            ChangeEvent cev = new ChangeEvent(this, Annotation.PROPERTY, key);
            featureAnnotationHub.firePreChange(cev);
            _setProperty(key, value);
            featureAnnotationHub.firePostChange(cev);
        }
    }

    private void _setProperty(Object key, Object value) 
        throws ChangeVetoException
    {
        if (underlyingAnnotation != null) {
            underlyingAnnotation.setProperty(key, value);
        }
        persistProperty(key, value);
    }
    
    public void removeProperty(Object key)
        throws ChangeVetoException
    {
        if (underlyingAnnotation == null) {
            initAnnotations();
        }
        if (!underlyingAnnotation.containsProperty(key)) {
            throw new NoSuchElementException("Annotation doesn't contain property " + key.toString());
        }
        
        BioSQLFeatureAnnotationChangeHub featureAnnotationHub = seqDB.getFeatureAnnotationChangeHub();
        synchronized (featureAnnotationHub) {
            ChangeEvent cev = new ChangeEvent(this, Annotation.PROPERTY, key);
            featureAnnotationHub.firePreChange(cev);
            underlyingAnnotation.removeProperty(key);
            persistProperty(key, null);
            featureAnnotationHub.firePostChange(cev);
        }
    }

    private void persistProperty(Object key, Object value)
        throws ChangeVetoException
    {
	Connection conn = null;
	try {
	    conn = seqDB.getDataSource().getConnection();
	    conn.setAutoCommit(false);

	    seqDB.getFeaturesSQL().persistProperty(conn, feature_id, key, value, true);

	    conn.commit();
	    conn.close();
	} catch (SQLException ex) {
	    boolean rolledback = false;
	    if (conn != null) {
		try {
		    conn.rollback();
		    rolledback = true;
		} catch (SQLException ex2) {}
            try {conn.close();} catch (SQLException ex3) {}
	    }
	    throw new BioRuntimeException("Error adding BioSQL tables" + 
					(rolledback ? " (rolled back successfully)" : ""), ex);
	}
    }

    public boolean containsProperty(Object key) {
	if (underlyingAnnotation == null) {
	    initAnnotations();
	}

	return underlyingAnnotation.containsProperty(key);
    }

    public Set keys() {
	if (underlyingAnnotation == null) {
	    initAnnotations();
	}

	return underlyingAnnotation.keys();
    }

    public Map asMap() {
	if (underlyingAnnotation == null) {
	    initAnnotations();
	}

	return Collections.unmodifiableMap(underlyingAnnotation.asMap());
    }
    
    public void addChangeListener(ChangeListener cl) {
	addChangeListener(cl, ChangeType.UNKNOWN);
    }
    
    public void addChangeListener(ChangeListener cl, ChangeType ct) {
	seqDB.getFeatureAnnotationChangeHub().addListener(new Integer(feature_id), cl, ct);
    }

    public void removeChangeListener(ChangeListener cl) {
	removeChangeListener(cl, ChangeType.UNKNOWN);
    }

    public void removeChangeListener(ChangeListener cl, ChangeType ct) {
	seqDB.getFeatureAnnotationChangeHub().removeListener(new Integer(feature_id), cl, ct);
    }

    public boolean isUnchanging(ChangeType ct) {
	return false;
    }
}
