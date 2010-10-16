/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.symbol.FuzzyLocation;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;

/**
 * Behind-the-scenes adaptor to the features sub-schema of BioSQL.
 *
 * @author Thomas Down
 * @author Simon Foote
 * @author Len Trigg
 * @author Richard Holland
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */
class FeaturesSQL {
    private BioSQLSequenceDB seqDB;
		private HashMap rankType;

    FeaturesSQL(BioSQLSequenceDB seqDB) {
	this.seqDB = seqDB;
    }
    
    //
    // Feature retrieval
    //

    /**
     * Get some features out of BioSQL and fire SeqIO events to the specified listener.
     * Currently 4 modes:
     *
     * <ul>
     * <li>Get all features on a bioentry (including children)</li>
     * <li>Get all top-level features in a region</li>
     * <li>Get all children of a specified parent</li>
     * <li>Get a particular feature by ID</li>
     * </ul>
     *
     * This is rather ugly, but it's well hidden.  Not sure what the options would be for
     * a cleaner API.
     */

    public void retrieveFeatures(int bioentry_id, 
				                 SeqIOListener listener,
                                 Location overlappingRegion,
                                 int immediateChildrenOfParent,
                                 int featureID) 
        throws SQLException, BioException
    {
        Connection conn = seqDB.getDataSource().getConnection();
        Map fmap = new HashMap();
        Map qmap = new HashMap();
        Map lmap = new HashMap();

	
        PreparedStatement get_features = null;
        if (overlappingRegion == null && immediateChildrenOfParent < 0 && featureID < 0) {
            get_features = conn.prepareStatement(
			        "select seqfeature.seqfeature_id, " +
                    "       seqfeature.type_term_id, " +
                    "       seqfeature.source_term_id " +
                    "  from seqfeature " +
                    " where seqfeature.bioentry_id = ?"
			);
            get_features.setInt(1, bioentry_id);
        } else if (overlappingRegion != null) {
            get_features = conn.prepareStatement(
			        "select seqfeature.seqfeature_id, " +
                    "       seqfeature.type_term_id, " +
                    "       seqfeature.source_term_id " +
                    "  from location, seqfeature " +
                    " where seqfeature.bioentry_id = ? and " +
                    "       location.seqfeature_id = seqfeature.seqfeature_id and " +
                    "       location.end_pos >= ? and " +
                    "       location.start_pos <= ? " +
                    " group by seqfeature.seqfeature_id, seqfeature.type_term_id, seqfeature.source_term_id"
             );
             get_features.setInt(1, bioentry_id);
             get_features.setInt(2, overlappingRegion.getMin());
             get_features.setInt(3, overlappingRegion.getMax());
        } else if (immediateChildrenOfParent >= 0) {
            get_features = conn.prepareStatement(
			        "select seqfeature.seqfeature_id, " +
                    "       seqfeature.type_term_id, " +
                    "       seqfeature.source_term_id " +
                    "  from seqfeature, seqfeature_relationship " +
                    " where seqfeature.seqfeature_id = seqfeature_relationship.subject_seqfeature_id and " +
                    "       seqfeature_relationship.object_seqfeature_id = ?"
	         );
             get_features.setInt(1, immediateChildrenOfParent);
        } else if (featureID >= 0) {
	        get_features = conn.prepareStatement(
			        "select seqfeature.seqfeature_id, " +
                    "       seqfeature.type_term_id, " +
                    "       seqfeature.source_term_id, " +
                    "       seqfeature.bioentry_id " + 
                    "  from seqfeature " +
                    " where seqfeature.seqfeature_id = ?"
	        );
            get_features.setInt(1, featureID);
        } else {
            conn.close();
            throw new BioException("I'm afraid you can't do that!");
        }

        ResultSet rs = get_features.executeQuery();
        while (rs.next()) {
            int feature_id = rs.getInt(1);
            StrandedFeature.Template templ = new StrandedFeature.Template();
            templ.type = seqDB.getOntologyTerm(rs.getInt(2));
            templ.source = seqDB.getOntologyTerm(rs.getInt(3));
            templ.annotation = new BioSQLFeatureAnnotation(seqDB, feature_id);
            // templ.annotation = new SmallAnnotation();
            fmap.put(new Integer(feature_id), templ);

            if (featureID >= 0 && bioentry_id < 0) {
                bioentry_id = rs.getInt(4);
                listener.addSequenceProperty("_biosql_internal.bioentry_id", new Integer(bioentry_id));
            }
        }
        rs.close();
        get_features.close();

        // Fetch annotations (worth a try!)
        
        PreparedStatement get_annotations = null;
        if (overlappingRegion == null && immediateChildrenOfParent < 0 && featureID < 0) {
            get_annotations = conn.prepareStatement(
		            "select seqfeature_qualifier_value.seqfeature_id, " +
                    "       seqfeature_qualifier_value.term_id, " +
                    "       seqfeature_qualifier_value.value " +
                    "  from seqfeature, seqfeature_qualifier_value " +
                    " where seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id and " +
                    "       seqfeature.bioentry_id = ?"
            );
            get_annotations.setInt(1, bioentry_id);
        } else if (overlappingRegion != null) {
           get_annotations = conn.prepareStatement(
		            "select seqfeature_qualifier_value.seqfeature_id, " +
                    "       seqfeature_qualifier_value.term_id, " +
                    "       seqfeature_qualifier_value.value " +
                    "  from seqfeature, seqfeature_qualifier_value, location " +
                    " where seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id and " +
                    "       seqfeature.bioentry_id = ? and" +
                    "       location.seqfeature_id = seqfeature.seqfeature_id and " +
                    "       location.end_pos >= ? and " +
                    "       location.start_pos <= ? " +
                    "       group by seqfeature_qualifier_value.seqfeature_id, " +
                    "                seqfeature_qualifier_value.term_id, " +
                    "                seqfeature_qualifier_value.value"
           );
           get_annotations.setInt(1, bioentry_id);
           get_annotations.setInt(2, overlappingRegion.getMin());
           get_annotations.setInt(3, overlappingRegion.getMax());
        } else if (immediateChildrenOfParent >= 0) {
            get_annotations = conn.prepareStatement(
           		   "select seqfeature_qualifier_value.seqfeature_id, " +
                   "       seqfeature_qualifier_value.term_id, " +
                   "       seqfeature_qualifier_value.value " +
                   "  from seqfeature_qualifier_value, seqfeature_relationship " +
                   " where seqfeature_qualifier_value.seqfeature_id = seqfeature_relationship.subject_seqfeature_id and " +
                   "       seqfeature_relationship.object_seqfeature_id = ?"
           );
           get_annotations.setInt(1, immediateChildrenOfParent);
        } else if (featureID >= 0) {
           get_annotations = conn.prepareStatement(
		            "select seqfeature_qualifier_value.seqfeature_id, " +
                    "       seqfeature_qualifier_value.term_id, " +
                    "       seqfeature_qualifier_value.value " +
                    "  from seqfeature_qualifier_value " +
                    " where seqfeature_qualifier_value.seqfeature_id = ?"
            );
           get_annotations.setInt(1, featureID);
        }
        rs = get_annotations.executeQuery();
        while (rs.next()) {
            Integer fid = new Integer(rs.getInt(1));
            String key = seqDB.getOntologyTerm(rs.getInt(2));
            String value = rs.getString(3).trim();
            Feature.Template templ = (Feature.Template) fmap.get(fid);
            try {
                ((BioSQLFeatureAnnotation) templ.annotation).initProperty(key, value);
            } catch (ChangeVetoException ex) {
	          try {conn.close();} catch (SQLException ex3) {}
                throw new BioError("Couldn't modify hidden FeatureHolder");
            }
        }
        rs.close();
        get_annotations.close();

	// Fetch those crappy location qualifiers first...

	/*

	if (seqDB.isLocationQualifierSupported()) {
	    PreparedStatement get_location_crap = conn.prepareStatement(
			    "select location_qualifier_value.location_id, " +
			    "       seqfeature_qualifier.qualifier_name, " +
			    "       location_qualifier_value.qualifier_value, " +
			    "       location_qualifier_value.qualifier_int_value " +
			    "  from location_qualifier_value, location, seqfeature, seqfeature_qualifier " +
			    " where seqfeature.bioentry_id = ? and " +
			    "       location.seqfeature_id = seqfeature.seqfeature_id and " +
			    "       location_qualifier_value.location_id = location.location_id and " +
			    "       seqfeature_qualifier.seqfeature_qualifier_id = location_qualifier_value.seqfeature_qualifier_id");
	    get_location_crap.setInt(1, bioentry_id);
	    rs = get_location_crap.executeQuery();
	    while (rs.next()) {
		LocationQualifierMemento lqm = new LocationQualifierMemento();
		int location_id = rs.getInt(1);
		lqm.qualifier_name = rs.getString(2).trim();    // HACK due to stupid schema change
		lqm.qualifier_value = rs.getString(3).trim();
		lqm.qualifier_int = rs.getInt(4);
		
		Integer location_id_boxed = new Integer(location_id);
		List l = (List) qmap.get(location_id_boxed);
		if (l == null) {
		    l = new ArrayList();
		    qmap.put(location_id_boxed, l);
		}
		l.add(lqm);
	    }
	}

	*/

	
	    // Fetch locations
       
        PreparedStatement get_locations = null;
        if (overlappingRegion == null && immediateChildrenOfParent < 0 && featureID < 0) {
            get_locations = conn.prepareStatement(
		            "select location.location_id, " +
                    "       location.seqfeature_id, " +
                    "       location.start_pos, " +
                    "       location.end_pos, " +
                    "       location.strand " +
                    "  from seqfeature, location " +
                    " where location.seqfeature_id = seqfeature.seqfeature_id and " +
                    "       seqfeature.bioentry_id = ?"
            );
            get_locations.setInt(1, bioentry_id);
        } else if (overlappingRegion != null) {
            get_locations = conn.prepareStatement(
		           "select location.location_id, " +
                   "       location.seqfeature_id, " +
                   "       location.start_pos, " +
                   "       location.end_pos, " +
                   "       location.strand " +
                   "  from location, location as sfl2, seqfeature " +
                   " where location.seqfeature_id = seqfeature.seqfeature_id and " +
                   "       seqfeature.bioentry_id = ? and " +
                   "       sfl2.seqfeature_id = seqfeature.seqfeature_id and " +
                   "       sfl2.end_pos >= ? and " +
                   "       sfl2.start_pos <= ? " +
                   " group by location.location_id, " +
                   "          location.seqfeature_id, " +
                   "          location.start_pos, " +
                   "          location.end_pos, " +
                   "          location.strand"
            );
            get_locations.setInt(1, bioentry_id);
            get_locations.setInt(2, overlappingRegion.getMin());
            get_locations.setInt(3, overlappingRegion.getMax());
        } else if (immediateChildrenOfParent >= 0) {
            get_locations = conn.prepareStatement(
		           "select location.location_id, " +
                   "       location.seqfeature_id, " +
                   "       location.start_pos, " +
                   "       location.end_pos, " +
                   "       location.strand " +
                   "  from location, seqfeature_relationship " +
                   " where location.seqfeature_id = seqfeature_relationship.subject_seqfeature_id and " +
                   "       seqfeature_relationship.object_seqfeature_id = ?"
            );
            get_locations.setInt(1, immediateChildrenOfParent);
        } else if (featureID >= 0) {
            get_locations = conn.prepareStatement(
		            "select location.location_id, " +
                    "       location.seqfeature_id, " +
                    "       location.start_pos, " +
                    "       location.end_pos, " +
                    "       location.strand " +
                    "  from location " +
                    " where location.seqfeature_id = ?");
            get_locations.setInt(1, featureID);
        }


        rs = get_locations.executeQuery();
        while (rs.next()) {
            Integer lid = new Integer(rs.getInt(1));
            Integer fid = new Integer(rs.getInt(2));
            int start = rs.getInt(3);
            int end = rs.getInt(4);
            int istrand = rs.getInt(5);
            
            StrandedFeature.Strand strand = StrandedFeature.UNKNOWN;
            if (istrand > 0) {
                strand = StrandedFeature.POSITIVE;
            } else if (istrand < 0) {
                strand = StrandedFeature.NEGATIVE;
            }
            StrandedFeature.Template templ = (StrandedFeature.Template) fmap.get(fid);
            if (templ.strand != null && templ.strand != strand) {
                // throw new BioRuntimeException("Feature strands don't match");
                // Really don't want to support these at all, but...
                templ.strand = StrandedFeature.UNKNOWN;
            } else {
                templ.strand = strand;
            }
           
            Location bloc;
            if (start == end) {
                bloc = new PointLocation(start);
            } else {
                bloc = new RangeLocation(start, end);
            }
           
            List locationCrap = (List) qmap.get(lid);
            if (locationCrap != null) {
                int min_start = -1;
                int min_end = -1;
                int max_start = -1;
                int max_end = -1;
                boolean unknown_start = false;
                boolean unknown_end = false;
                boolean unbounded_start = false;
                boolean unbounded_end = false;
                boolean isFuzzy = false;
		
                for (Iterator i = locationCrap.iterator(); i.hasNext(); ) {
                    LocationQualifierMemento lqm = (LocationQualifierMemento) i.next();
                    String qname = lqm.qualifier_name;
		    
                    if ("min_start".equals(qname)) {
                        min_start = lqm.qualifier_int;
                        isFuzzy = true;
                    } else if ("max_start".equals(qname)) {
                        max_start = lqm.qualifier_int;
                        isFuzzy = true;
                    } else if ("min_end".equals(qname)) {
                        min_end = lqm.qualifier_int;
                        isFuzzy = true;
                    } else if ("max_end".equals(qname)) {
                        max_end = lqm.qualifier_int;
                        isFuzzy = true;
                    } else if ("start_pos_type".equals(qname)) {
                        if ("BEFORE".equalsIgnoreCase(lqm.qualifier_value)) {
                            unbounded_start = true;
                            isFuzzy = true;
                        }
                    } if ("end_pos_type".equals(qname)) {
                        if ("AFTER".equalsIgnoreCase(lqm.qualifier_value)) {
                            unbounded_end = true;
                            isFuzzy = true;
                        }
                    } 
                }

                if (isFuzzy) {
                    if (unknown_start) {
                        min_start = Integer.MIN_VALUE;
                        max_start = Integer.MAX_VALUE;
                    }
                    if (unbounded_start) {
                        min_start = Integer.MIN_VALUE;
                    }
                    if (unknown_end) {
                        min_end = Integer.MIN_VALUE;
                        max_end = Integer.MAX_VALUE;
                    }
                    if (unbounded_end) {
                        max_end = Integer.MAX_VALUE;
                    }
		    
                    if (min_start == -1) {
                        min_start = bloc.getMin();
                    }
                    if (max_start == -1) {
                        max_start = bloc.getMin();
                    }
                    if (min_end == -1) {
                        min_end = bloc.getMax();
                    } 
                    if (max_end == -1) {
                        max_end = bloc.getMax();
                    }
		    
                    bloc = new FuzzyLocation(min_start,
                                             max_end,
                                             max_start,
                                             min_end,
                                             FuzzyLocation.RESOLVE_INNER);
                }
            }

            List ll = (List) lmap.get(fid);
            if (ll == null) {
                ll = new ArrayList();
                lmap.put(fid, ll);
            }
            ll.add(bloc);
        }
        rs.close();
        get_locations.close();
	
        // Bind location information to features
	
    	for (Iterator i = fmap.entrySet().iterator(); i.hasNext(); ) {
            Map.Entry me = (Map.Entry) i.next();
            Integer fid = (Integer) me.getKey();
            StrandedFeature.Template templ = (StrandedFeature.Template) me.getValue();
	    
    	    List ll = (List) lmap.get(fid);
            if (ll == null) {
                templ.location = Location.empty;
                
	        //  conn.close();
                //throw new BioRuntimeException("BioSQL SeqFeature doesn't have any associated location spans. seqfeature_id=" + fid);

            } else {
	    
    	      Location loc = null;
              if (ll.size() == 1) {
                  loc = (Location) ll.get(0);
              } else {
                  loc = LocationTools.union(ll);
              }
              templ.location = loc;
            }
        }

        // Check hierarchy
	
    	Set toplevelFeatures = new HashSet(fmap.keySet());
        Map featureHierarchy = new HashMap();
        int specifiedParent = -1;
        if (immediateChildrenOfParent < 0 && featureID < 0) {
            PreparedStatement get_hierarchy;
            if (overlappingRegion == null) {
                get_hierarchy = conn.prepareStatement(
                        "select object_seqfeature_id, subject_seqfeature_id " +
                        "  from seqfeature_relationship, seqfeature " +
                        " where object_seqfeature_id = seqfeature.seqfeature_id and " +
                        "       seqfeature.bioentry_id = ?"
                );
                get_hierarchy.setInt(1, bioentry_id);
            } else {
                get_hierarchy = conn.prepareStatement(
                        "select object_seqfeature_id, subject_seqfeature_id " +
                        "  from seqfeature_relationship, seqfeature, location " +
                        " where object_seqfeature_id = seqfeature.seqfeature_id and " +
                        "       seqfeature.bioentry_id = ? and " +
                        "       location.seqfeature_id = object_seqfeature_id and " +
                        "       location.end_pos >= ? and " +
                        "       location.start_pos <= ? " +
                        "       group by object_seqfeature_id, subject_seqfeature_id"
                 );
                 get_hierarchy.setInt(1, bioentry_id);
                 get_hierarchy.setInt(2, overlappingRegion.getMin());
                 get_hierarchy.setInt(3, overlappingRegion.getMax());
            }
            rs = get_hierarchy.executeQuery();
            while (rs.next()) {
                Integer parent = new Integer(rs.getInt(1));
                Integer child = new Integer(rs.getInt(2));
		
                toplevelFeatures.remove(child);
                List cl = (List) featureHierarchy.get(parent);
                if (cl == null) {
                    cl = new ArrayList();
                    featureHierarchy.put(parent, cl);
                }
                cl.add(child);
            }
            rs.close();
            get_hierarchy.close();
        } else if (immediateChildrenOfParent >= 0) {
            specifiedParent = immediateChildrenOfParent;
        } else if (featureID >= 0) {
            PreparedStatement discover_parent = conn.prepareStatement(
		            "select object_seqfeature_id " +
                    "  from seqfeature_relationship " +
                    " where object_seqfeature_id = ?"
            );
            discover_parent.setInt(1, featureID);
            rs = discover_parent.executeQuery();
            if (rs.next()) {
                specifiedParent = rs.getInt(1);
            }
            rs.close();
            discover_parent.close();
        }

        conn.close();
        conn = null;

        for (Iterator tlfi = toplevelFeatures.iterator(); tlfi.hasNext(); ) {
            Integer fid = (Integer) tlfi.next();
            Feature.Template templ = (Feature.Template) fmap.get(fid);
            boolean childrenFetched = (immediateChildrenOfParent < 0);
            if (overlappingRegion != null) {
                if (!overlappingRegion.contains(templ.location)) {
                    childrenFetched = false;
                }
            }
            fireFeatureTree(listener, 
                            fid,
                            fmap,
                            featureHierarchy,
                            childrenFetched,
                            new Integer(specifiedParent)
                           );
        }
    }

    private void fireFeatureTree(SeqIOListener listener,
				 Integer fid,
				 Map fmap,
				 Map featureHierarchy,
				 boolean childrenFetched,
				 Integer pid)
        throws BioException
    {
        Feature.Template templ = (Feature.Template) fmap.get(fid);
        listener.startFeature(templ);
        listener.addFeatureProperty("_biosql_internal.feature_id", fid);
        listener.addFeatureProperty("_biosql_internal.parent_id", pid);
        if (childrenFetched) {
            List children = (List) featureHierarchy.get(fid);
            if (children == null) {
                listener.addFeatureProperty("_biosql_internal.hint_childfree", Boolean.TRUE);
            } else {
                for (Iterator ci = children.iterator(); ci.hasNext(); ) {
                    Integer childID = (Integer) ci.next();
                    fireFeatureTree(listener, childID, fmap, featureHierarchy, childrenFetched, fid);
                }
            }
        }
        listener.endFeature();
    }

    
    private static class LocationQualifierMemento {
	public String qualifier_name;
	public String qualifier_value;
	public int qualifier_int;
    }

    //
    // Feature live updates
    //

    void setFeatureType(int feature_id, String type)
        throws SQLException
    {
	Connection conn = null;
	try {
	    conn = seqDB.getDataSource().getConnection();
	    conn.setAutoCommit(false);
	    
	    int seqfeature_key = seqDB.intern_ontology_term(conn, type);
	    PreparedStatement update_key = conn.prepareStatement("update seqfeature " + 
                                                                 "   set type_term_id = ? " +
								 " where seqfeature_id = ?");
	    update_key.setInt(1, seqfeature_key);
	    update_key.setInt(2, feature_id);
	    update_key.executeUpdate();
            update_key.close();
	    conn.commit();
	    conn.close();
	} catch (SQLException ex) {
	    if (conn != null) {
		try {
		    conn.rollback();
		} catch (SQLException ex2) {}
	          try {conn.close();} catch (SQLException ex3) {}
	    }
	    throw ex;
	}
    }

    void setFeatureSource(int feature_id, String source)
        throws SQLException
    {
	Connection conn = null;
	try {
	    conn = seqDB.getDataSource().getConnection();
	    conn.setAutoCommit(false);
	    
	    int seqfeature_source = seqDB.intern_ontology_term(conn, source);
	    PreparedStatement update_source = conn.prepareStatement("update seqfeature " + 
								    "   set source_term_id = ? " +
								    " where seqfeature_id = ?");
	    update_source.setInt(1, seqfeature_source);
	    update_source.setInt(2, feature_id);
	    update_source.executeUpdate();
            update_source.close();
	    conn.commit();
	    conn.close();
	} catch (SQLException ex) {
	    if (conn != null) {
		try {
		    conn.rollback();
		} catch (SQLException ex2) {}
                try {conn.close();} catch (SQLException ex3) {}
	    }
	    throw ex;
	}
    }

    void setFeatureLocation(int feature_id, Location location, StrandedFeature.Strand s)
        throws SQLException
    {
	Connection conn = null;
	try {
	    conn = seqDB.getDataSource().getConnection();
	    conn.setAutoCommit(false);
	    
	    PreparedStatement del_oldlocation = conn.prepareStatement(
		    "delete from location " +
		    " where seqfeature_id = ?"
	    );
	    del_oldlocation.setInt(1, feature_id);
	    del_oldlocation.executeUpdate();
	    del_oldlocation.close();

	    PreparedStatement add_locationspan = conn.prepareStatement(
                    "insert into location " +
	            "       (seqfeature_id, start_pos, end_pos, strand, rank) " +
    		    "values (?, ?, ?, ?, ?)"
	    );

	    int strandNum;
	    if (s == StrandedFeature.POSITIVE) {
		strandNum = 1;
	    } else if (s== StrandedFeature.NEGATIVE) {
		strandNum = -1;
	    } else {
		strandNum = 0;
	    }

	    int rank = 0;
	    for (Iterator i = location.blockIterator(); i.hasNext(); ) {
		Location bloc = (Location) i.next();
		add_locationspan.setInt(1, feature_id);
		add_locationspan.setInt(2, bloc.getMin());
		add_locationspan.setInt(3, bloc.getMax());
		add_locationspan.setInt(4, strandNum);
		add_locationspan.setInt(5, ++rank);
		add_locationspan.executeUpdate();
	    }
	    add_locationspan.close();

	    conn.commit();
	    conn.close();
	} catch (SQLException ex) {
	    if (conn != null) {
		try {
		    conn.rollback();
		} catch (SQLException ex2) {}
            try {conn.close();} catch (SQLException ex3) {}
	    }
	    throw ex;
	}
    }

    //
    // Feature persistance
    //
    
    void persistFeatures(Connection conn, int bioentry_id, FeatureHolder features)
        throws BioException, SQLException
    {
        persistFeatures(conn, bioentry_id, features, -1);
    }

    private void persistFeatures(Connection conn, int bioentry_id, FeatureHolder features, int parent)
        throws BioException, SQLException
    {
			
	// HashMap to track rank for each feature_type
	// This allows the unique key constraint in seqfeature to be valid as
	// bioentries can have multiple features of same type/source (ie. Genbank CDS)
	// This way is faster than just increasing the rank value each time
	// Only set Map for parent, otherwise adding children of same type throws 
	// a duplicate key error
	if (parent < 0) {
		rankType = new HashMap();
	}
	
	for (Iterator fi = features.features(); fi.hasNext(); ) {
	    Feature f = (Feature) fi.next();
			
            // Get next rank value for feature type
            int rank = 0;
            String fType = f.getType();
            if (rankType.containsKey(fType)) {
                rank = ((Integer) rankType.get(fType)).intValue() + 1;
                rankType.put(fType, new Integer(rank));
            } else {
                rankType.put(fType, new Integer(0));
            }
            
	    if (! (f instanceof ComponentFeature)) {
		int id = persistFeature(conn, bioentry_id, f, parent, rank);
		if (seqDB.isHierarchySupported()) {
		    persistFeatures(conn, bioentry_id, f, id);
		}
	    }
	}

    }
		

    int persistFeature(Connection conn,
		       int bioentry_id,
		       Feature f,
		       int parent_id,
           int typeRank)
	throws BioException, SQLException
    {
	int id = -1;
	boolean locationWritten = false;

	if (seqDB.isSPASupported()) {
	    if (f.getLocation().isContiguous()) {
		Location loc = f.getLocation();

		PreparedStatement add_feature = conn.prepareStatement(
		        "select create_seqfeature_onespan(?, ?, ?, ?, ?, ?)"
		);
		add_feature.setInt(1, bioentry_id);
		add_feature.setString(2, f.getType());
		add_feature.setString(3, f.getSource());
		add_feature.setInt(4, loc.getMin());
		add_feature.setInt(5, loc.getMax());
		if (f instanceof StrandedFeature) {
		    StrandedFeature.Strand s = ((StrandedFeature) f).getStrand();
		    if (s == StrandedFeature.POSITIVE) {
			add_feature.setInt(6, 1);
		    } else if (s== StrandedFeature.NEGATIVE) {
			add_feature.setInt(6, -1);
		    } else {
			add_feature.setInt(6, 0);
		    }
		} else {
		    add_feature.setInt(6, 0);
		}
		ResultSet rs = add_feature.executeQuery();
		if (rs.next()) {
		    id = rs.getInt(1);
		}
                rs.close();
		add_feature.close();

		locationWritten = true;
	    } else {
		PreparedStatement add_feature = conn.prepareStatement(
		        "select create_seqfeature(?, ?, ?)"
		);
		add_feature.setInt(1, bioentry_id);
		add_feature.setString(2, f.getType());
		add_feature.setString(3, f.getSource());
		ResultSet rs = add_feature.executeQuery();
		if (rs.next()) {
		    id = rs.getInt(1);
		}
                rs.close();
		add_feature.close();
	    }
	} else {
	    int seqfeature_key = seqDB.intern_ontology_term(conn, f.getType());
	    int seqfeature_source = seqDB.intern_ontology_term(conn, f.getSource());
			// Because of unique key constraints on seqfeature
			// Need to select the maximum rank value for the bioentry,key,source value
			// if rank is set to -1
			if (typeRank < 0) {
				PreparedStatement select_rank = conn.prepareStatement(
					"select max(rank) from seqfeature where bioentry_id=?"
					+ " and type_term_id=? and source_term_id=?");
				select_rank.setInt(1, bioentry_id);
				select_rank.setInt(2, seqfeature_key);
				select_rank.setInt(3, seqfeature_source);
				ResultSet rs = select_rank.executeQuery();
				if (rs.next()) {
					typeRank = rs.getInt(1) + 1;
				}
                                rs.close();
                                select_rank.close();
			}
			
	    PreparedStatement add_feature = conn.prepareStatement(
		"insert into seqfeature "+
		"       (bioentry_id, type_term_id, source_term_id, rank) " +
		"values (?, ?, ?, ?)"
	    );
	    add_feature.setInt(1, bioentry_id);
	    add_feature.setInt(2, seqfeature_key);
	    add_feature.setInt(3, seqfeature_source);
            add_feature.setInt(4, typeRank);
	    add_feature.executeUpdate();
	    add_feature.close();

	    id = seqDB.getDBHelper().getInsertID(conn, "seqfeature", "seqfeature_id");
	}

	if (!locationWritten) {
	    PreparedStatement add_locationspan = conn.prepareStatement(
                    "insert into location " +
	            "       (seqfeature_id, start_pos, end_pos, strand, rank) " +
    		    "values (?, ?, ?, ?, ?)"
	    );

	    int strandNum;

	    if (f instanceof StrandedFeature) {
		StrandedFeature.Strand s = ((StrandedFeature) f).getStrand();
		if (s == StrandedFeature.POSITIVE) {
		    strandNum = 1;
		} else if (s== StrandedFeature.NEGATIVE) {
		    strandNum = -1;
		} else {
		    strandNum = 0;
		}
	    } else {
		strandNum = 0;
	    }

	    int rank = 0;
	    for (Iterator i = f.getLocation().blockIterator(); i.hasNext(); ) {
		Location bloc = (Location) i.next();
		add_locationspan.setInt(1, id);
		add_locationspan.setInt(2, bloc.getMin());
		add_locationspan.setInt(3, bloc.getMax());
		add_locationspan.setInt(4, strandNum);
		add_locationspan.setInt(5, ++rank);
		add_locationspan.executeUpdate();
	    }
	    add_locationspan.close();
	}

	//
	// Persist anything in the annotation bundle, as well.
	//

	for (Iterator ai = f.getAnnotation().asMap().entrySet().iterator(); ai.hasNext(); ) {
	    Map.Entry akv = (Map.Entry) ai.next();
	    persistProperty(conn, id, akv.getKey(), akv.getValue(), false);
	}

	//
	// Persist link to parent
	//

	if (parent_id >= 0) {
	    PreparedStatement add_hierarchy = conn.prepareStatement(
		"insert into seqfeature_relationship "+
		"       (object_seqfeature_id, subject_seqfeature_id, term_id) " +
		"values (?, ?, ?)"
		);
	    add_hierarchy.setInt(1, parent_id);
	    add_hierarchy.setInt(2, id);
	    add_hierarchy.setInt(3, seqDB.intern_ontology_term(conn, "contains"));
	    add_hierarchy.executeUpdate();
	    add_hierarchy.close();
	}

	return id;
    }

    void removeFeature(BioSQLFeature f)
        throws ChangeVetoException
    {
        Connection conn = null;
        try {
            conn = seqDB.getDataSource().getConnection();
            conn.setAutoCommit(false);

            removeFeature(conn, f);

            conn.commit();
            conn.close();
        } catch (SQLException ex) {
	    boolean rolledback = false;
	    if (conn != null) {
		try {
		    conn.rollback();
		    rolledback = true;
		} catch (SQLException ex2) {}
	    }
	    throw new BioRuntimeException("Error removing from BioSQL tables" 
					+ (rolledback ? " (rolled back successfully)" : ""), ex);
	}
    }

    private void removeFeature(Connection conn, BioSQLFeature f)
        throws SQLException, ChangeVetoException
    {
        Iterator children = ((FeatureHolder) f).features();
        while (children.hasNext()) {
            Feature f2 = (Feature) children.next();
            if (f2 instanceof BioSQLFeature) {
                removeFeature(conn, (BioSQLFeature) f2);
            }
        }

        int feature_id = f._getInternalID();

        PreparedStatement delete_locs = conn.prepareStatement("delete from location " +
                                                              " where location.seqfeature_id = ?");
        delete_locs.setInt(1, feature_id);
        delete_locs.executeUpdate();
        delete_locs.close();

        PreparedStatement delete_fqv = conn.prepareStatement("delete from seqfeature_qualifier_value " +
                                                             " where seqfeature_qualifier_value.seqfeature_id = ?");
        delete_fqv.setInt(1, feature_id);
        delete_fqv.executeUpdate();
        delete_fqv.close();

        PreparedStatement delete_rel = conn.prepareStatement("delete from seqfeature_relationship " +
                                                             " where subject_seqfeature_id = ?");
        delete_rel.setInt(1, feature_id);
        delete_rel.executeUpdate();
        delete_rel.close();

        PreparedStatement delete_feature = conn.prepareStatement("delete from seqfeature " +
                                                                 " where seqfeature_id = ?");
        delete_feature.setInt(1, feature_id);
        delete_feature.executeUpdate();
        delete_feature.close();
    }

    /**
     * Persist a property.  Nothing is written if value is void
     */
    
    void persistProperty(Connection conn,
			 int feature_id,
			 Object key,
			 Object value,
			 boolean removeFirst)
        throws SQLException
    {
	String keyString = key.toString();

	if (removeFirst) {
	    int id = seqDB.intern_ontology_term(conn, keyString);
	    PreparedStatement remove_old_value = conn.prepareStatement("delete from seqfeature_qualifier_value " +
								       " where seqfeature_id = ? and term_id = ?");
	    remove_old_value.setInt(1, feature_id);
	    remove_old_value.setInt(2, id);
	    remove_old_value.executeUpdate();
	    remove_old_value.close();
	}

        if (value != null) {
            PreparedStatement insert_new;
            if (seqDB.isSPASupported()) {
                insert_new = conn.prepareStatement("insert into seqfeature_qualifier_value " +
                                                   "       (seqfeature_id, term_id, rank, value) " +
                                                   "values (?, intern_ontology_term( ? ), ?, ?)");
                if (value instanceof Collection) {
                    int cnt = 0;
                    for (Iterator i = ((Collection) value).iterator(); i.hasNext(); ) {
                        insert_new.setInt(1, feature_id);
                        insert_new.setString(2, keyString);
                        insert_new.setInt(3, ++cnt);
                        insert_new.setString(4, i.next().toString());
                        insert_new.executeUpdate();
                    }
                } else {
                    insert_new.setInt(1, feature_id);
                    insert_new.setString(2, keyString);
                    insert_new.setInt(3, 1);
                    insert_new.setString(4, value.toString());
                    insert_new.executeUpdate();
                }
                insert_new.close();
            } else {
                insert_new = conn.prepareStatement("insert into seqfeature_qualifier_value " +
                                                   "       (seqfeature_id, term_id, rank, value) " +
                                                   "values (?, ?, ?, ?)");
	        int sfq = seqDB.intern_ontology_term(conn, keyString);
                if (value instanceof Collection) {
                    int cnt = 0;
                    for (Iterator i = ((Collection) value).iterator(); i.hasNext(); ) {
                        insert_new.setInt(1, feature_id);
                        insert_new.setInt(2, sfq);
                        insert_new.setInt(3, ++cnt);
                        insert_new.setString(4, i.next().toString());
                        insert_new.executeUpdate();
                    }
                } else {
                    insert_new.setInt(1, feature_id);
                    insert_new.setInt(2, sfq);
                    insert_new.setInt(3, 1);
                    insert_new.setString(4, value.toString());
                    insert_new.executeUpdate();
                }
                insert_new.close();
            }
        }

    }

    void persistFeature(Feature f, int parent_id, int bioentry_id)
        throws BioException
    {
	Connection conn = null;
	try {
	    conn = seqDB.getDataSource().getConnection();
	    conn.setAutoCommit(false);
            // Set rank to -1, so will get looked up before feature added
	    int f_id = seqDB.getFeaturesSQL().persistFeature(conn, bioentry_id, f, parent_id, -1);
	    if (f instanceof BioSQLFeature) {
		((BioSQLFeature) f)._setInternalID(f_id);
		((BioSQLFeature) f)._setAnnotation(new BioSQLFeatureAnnotation(seqDB, f_id));
	    }
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
	    throw new BioException("Error adding to BioSQL tables" 
					+ (rolledback ? " (rolled back successfully)" : ""), ex);
	}
    }
}
