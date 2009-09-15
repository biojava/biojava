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

package org.biojavax.bio.seq;

import java.util.Set;
import java.util.TreeSet;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.ontology.InvalidTermException;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.RankedCrossRefable;
import org.biojavax.RichAnnotatable;
import org.biojavax.SimpleRichAnnotation;

/**
 * Represents a feature that can be given name and rank and be
 * moved from one sequence to another.
 * @author Richard Holland
 * @since 1.5
 */
public interface RichFeature extends StrandedFeature,RankedCrossRefable,RichAnnotatable,RichFeatureRelationshipHolder,Comparable {
    
    public static final ChangeType NAME = new ChangeType(
            "This feature's name has changed",
            "org.biojavax.bio.seq.RichFeature",
            "NAME"
            );
    public static final ChangeType RANK = new ChangeType(
            "This feature's rank has changed",
            "org.biojavax.bio.seq.RichFeature",
            "RANK"
            );
    public static final ChangeType SOURCETERM = new ChangeType(
            "This feature's source term has changed",
            "org.biojavax.bio.seq.RichFeature",
            "SOURCETERM"
            );
    public static final ChangeType TYPETERM = new ChangeType(
            "This feature's type term has changed",
            "org.biojavax.bio.seq.RichFeature",
            "TYPETERM"
            );
    public static final ChangeType LOCATION = new ChangeType(
            "This feature's location has changed",
            "org.biojavax.bio.seq.RichFeature",
            "LOCATION"
            );
    public static final ChangeType PARENT = new ChangeType(
            "This feature's parent has changed",
            "org.biojavax.bio.seq.RichFeature",
            "PARENT"
            );
    public static final ChangeType CROSSREF = new ChangeType(
            "This feature's crossrefs have changed",
            "org.biojavax.bio.seq.RichFeature",
            "CROSSREF"
            );
    public static final ChangeType RELATION = new ChangeType(
            "This feature's relations have changed",
            "org.biojavax.bio.seq.RichFeature",
            "RELATION"
            );
    
    /**
     * Sets the parent of this feature.
     * @param parent the parent the feature should identify itself with.
     * @throws ChangeVetoException if the new value is unacceptable.
     */
    public void setParent(FeatureHolder parent) throws ChangeVetoException;
    
    /**
     * Returns the name of this feature.
     * @return the name.
     */
    public String getName();
    
    /**
     * Sets the name of this feature.
     * @param name the name the feature should identify itself with.
     * @throws ChangeVetoException if the new value is unacceptable.
     */
    public void setName(String name) throws ChangeVetoException;
    
    /**
     * Returns the rank of this feature.
     * @return the rank.
     */
    public int getRank();
    
    /**
     * Sets the rank of this feature.
     * @param rank the rank the feature should identify itself with.
     * @throws ChangeVetoException if the new value is unacceptable.
     */
    public void setRank(int rank) throws ChangeVetoException;
    
    /**
     * Added-value extension of Feature.Template including bits we're 
     * interested in, eg. featureRelationshipSet for relationships with
     * other features, and rankedCrossRefs for lookups into other databases.
     */
    public static class Template extends Feature.Template {
        public Set featureRelationshipSet;
        public Set rankedCrossRefs;
    }
    
    /**
     * Some useful tools for working with features.
     */
    public static class Tools {
        
        // because it is a static class we don't want instances of it
        private Tools() {}
        
        /**
         * Takes a normal Feature and attempts to convert it into a RichFeature.
         * @param f the Feature to enrich.
         * @return the RichFeature equivalent.
         * @throws ChangeVetoException if any of the terms from the feature were
         * not convertible.
         */
        public static RichFeature enrich(Feature f) throws ChangeVetoException {
            try {
                if (f instanceof RichFeature) return (RichFeature)f;
                else return new SimpleRichFeature(f.getParent(),f.makeTemplate());
            } catch (InvalidTermException e) {
                throw new ChangeVetoException("Unable to convert one of the feature's terms",e);
            }
        }
        
        /**
         * A way to make a dummy feature
         */
        public static RichFeature makeEmptyFeature(){            
            RichFeature.Template templ = new RichFeature.Template();
            templ.annotation = new SimpleRichAnnotation();
            templ.featureRelationshipSet = new TreeSet();
            templ.rankedCrossRefs = new TreeSet();
            templ.type = "none";
            templ.source = "none";
            templ.location = RichLocation.full;
            try{
                return new SimpleRichFeature(new SimpleFeatureHolder(),templ);
            }catch(Exception ex){
                throw new BioError(ex);
            }//can't happen            
        }
    }
}

