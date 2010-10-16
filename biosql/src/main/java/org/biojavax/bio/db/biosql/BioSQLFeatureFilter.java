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

package org.biojavax.bio.db.biosql;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;

import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.utils.walker.WalkerFactory;
import org.biojavax.Note;
import org.biojavax.RichAnnotation;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleNote;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.RichLocation.Strand;
import org.biojavax.ontology.ComparableTerm;



/**
 * A filter for accepting or rejecting a feature.
 *
 * <p>
 * It is possible to write custom <code>FeatureFilter</code>s by implementing this
 * interface.  There are also a wide range of built-in features, and it is possible
 * to build complex queries using <code>FeatureFilter.And</code>, <code>FeatureFilter.Or</code>,
 * and <code>FeatureFilter.Not</code>.  Where possible, use of the built-in filters
 * is preferable to writing new filters, since the methods in the <code>FilterUtils</code>
 * class have access to special knowledge about the built-in filter types and how they
 * relate to one another.
 * </p>
 *
 * <p>
 * If the filter is to be used in a remote process, it is recognized that it may
 * be serialized and sent over to run remotely, rather than each feature being
 * retrieved locally.
 * </p>
 *
 * <p>
 * This class requires the Hibernate JAR files to be on your classpath at runtime. It is
 * designed ONLY for use with BioSQLRichSequenceDB and BioSQLBioEntryDB.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Richard Holland
 * @since 1.5
 * @since 1.5
 */

public interface BioSQLFeatureFilter extends FeatureFilter {
    
    /**
     * This method returns a Hibernate Criterion object that can be used to
     * query the database.
     * @return a Hibernate Criterion object representing this filter.
     */
    public Object asCriterion();
    
    /**
     * Returns a map of property names (keys) to aliases (values), if the criterion
     * returned by asCriterion() uses aliases at all. If not, then it must at least
     * return the empty map else you'll get NullPointerExceptions thrown elsewhere.
     * @return Map a map of property names to aliases used in the criterion.
     */
    public Map criterionAliasMap();
    
    /**
     * A class representing some useful stuff you can do with BioSQLFeatureFilters,
     * for instance converting plain FeatureFilters into a their BioSQLFeatureFilter
     * equivalents (where possible).
     */
    public static class Tools {
        /**
         * Convert a non-BioSQL FeatureFilter into a BioSQL one. We do this
         * by walking through it, converting any ones we recognise into their
         * BioSQLFeatureFilter equivalents. If we don't recognise them, we take
         * special action. For the child of an And, we can just ignore the missing
         * side and replace the And itself with the remaining side. For everything else,
         * the entire FeatureFilter is replaced by BioSQLFeatureFilter.all else we
         * run the risk of missing out potential candidates.
         * The end result is a filter that can be applied to the
         * database to filter out potential candidates for more rigorous selection
         * in-memory by the default filter() method in AbstractRichSequenceDB. Whether or
         * not the filter picks out everything correctly depends entirely on whether it
         * is made up of BioSQLFeatureFilter elements, or can be converted into them.
         */
        public static BioSQLFeatureFilter convert(FeatureFilter ff) {
            // The easy case first.
            if (ff instanceof BioSQLFeatureFilter) return (BioSQLFeatureFilter)ff;
            else {
                BioSQLFeatureFilter bff = attemptConversion(ff);
                if (bff!=null) return bff;
                else return BioSQLFeatureFilter.all; // catch-all case.
            }
        }
        
        private static BioSQLFeatureFilter attemptConversion(FeatureFilter ff) {
            // AND - convert both children. If both are convertible, return the And
            // of them. If only one is, return just that child. If neither are,
            // return null.
            if (ff instanceof FeatureFilter.And) {
                FeatureFilter.And ffand = (FeatureFilter.And)ff;
                BioSQLFeatureFilter child1 = attemptConversion(ffand.getChild1());
                BioSQLFeatureFilter child2 = attemptConversion(ffand.getChild2());
                if (child1==null && child2==null) return null;
                else if (child1==null && child2!=null) return child2;
                else if (child1!=null && child2==null) return child1;
                else return new BioSQLFeatureFilter.And(child1,child2);
            }
            // OR - convert both children. If both are convertible, return the Or
            // of them. Otherwise, return null.
            else if (ff instanceof FeatureFilter.Or) {
                FeatureFilter.Or ffor = (FeatureFilter.Or)ff;
                BioSQLFeatureFilter child1 = attemptConversion(ffor.getChild1());
                BioSQLFeatureFilter child2 = attemptConversion(ffor.getChild2());
                if (child1==null || child2==null) return null;
                else return new BioSQLFeatureFilter.Or(child1,child2);
            }
            // NOT - convert the child. If convertible, return the Not of it. Else,
            // return null.
            else if (ff instanceof FeatureFilter.Not) {
                FeatureFilter.Not ffnot = (FeatureFilter.Not)ff;
                BioSQLFeatureFilter child = attemptConversion(ffnot.getChild());
                if (child==null) return null;
                else return new BioSQLFeatureFilter.Not(child);
            }
            // BySource - convert the term to a Term from the default ontology then
            // try BySourceTerm.
            else if (ff instanceof FeatureFilter.BySource) {
                FeatureFilter.BySource ffsrc = (FeatureFilter.BySource)ff;
                String name = ffsrc.getSource();
                return new BioSQLFeatureFilter.BySourceTermName(name);
            }
            // ByType - convert the term to a Term from the default ontology then
            // try ByTypeTerm.
            else if (ff instanceof FeatureFilter.ByType) {
                FeatureFilter.ByType ffsrc = (FeatureFilter.ByType)ff;
                String name = ffsrc.getType();
                return new BioSQLFeatureFilter.ByTypeTermName(name);
            }
            // ContainedByLocation - simple pass-through
            else if (ff instanceof FeatureFilter.ContainedByLocation) {
                FeatureFilter.ContainedByLocation ffloc = (FeatureFilter.ContainedByLocation)ff;
                return new BioSQLFeatureFilter.ContainedByRichLocation(RichLocation.Tools.enrich(ffloc.getLocation()));
            }
            // BySequenceName - simple pass-through
            else if (ff instanceof FeatureFilter.BySequenceName) {
                FeatureFilter.BySequenceName ffsn = (FeatureFilter.BySequenceName)ff;
                return new BioSQLFeatureFilter.BySequenceName(ffsn.getSequenceName());
            }
            // ShadowOverlapsLocation - simple pass-through to OverlapsRichLocation, as we have no concept
            // of shadows within BioSQL so they are effectively the same thing.
            else if (ff instanceof FeatureFilter.ShadowOverlapsLocation) {
                FeatureFilter.ShadowOverlapsLocation ffloc = (FeatureFilter.ShadowOverlapsLocation)ff;
                return new BioSQLFeatureFilter.OverlapsRichLocation(RichLocation.Tools.enrich(ffloc.getLocation()));
            }
            // AnnotationContains - attempt to convert the key to a ComparableTerm, and the value to a string (retrieve the
            // sole member if it is a collection), then wrap the whole thing in a Note with rank 0 and try using
            // ByNoteWithValue to retrieve.
            else if (ff instanceof FeatureFilter.AnnotationContains) {
                FeatureFilter.AnnotationContains ffann = (FeatureFilter.AnnotationContains)ff;
                if (!(ffann.getValue() instanceof String)) return null;
                String noteValue = (String)ffann.getValue();
                ComparableTerm noteTerm;
                Object key = ffann.getKey();
                if (key instanceof Collection) {
                    Collection coll = (Collection)key;
                    if (coll.size()<1) return null;
                    else key = coll.toArray()[0];
                }
                if (key instanceof ComparableTerm) noteTerm = (ComparableTerm)key;
                else if (key instanceof String) noteTerm =  RichObjectFactory.getDefaultOntology().getOrCreateTerm((String)key);
                else return null;
                return new BioSQLFeatureFilter.ByNote(new SimpleNote(noteTerm,noteValue,0));
            }
            // StrandFilter - attempt to convert the StrandedFeature.Strand to a RichLocation.Strand then pass through
            // to ByStrand.
            else if (ff instanceof FeatureFilter.StrandFilter) {
                FeatureFilter.StrandFilter ffstr = (FeatureFilter.StrandFilter)ff;
                Strand strand = RichLocation.Strand.forName(""+ffstr.getStrand().getToken());
                return new BioSQLFeatureFilter.ByStrand(strand);
            }
            // HasAnnotation - attempt to convert the term into a ComparableTerm, then use ByNoteTermOnly.
            else if (ff instanceof FeatureFilter.HasAnnotation) {
                FeatureFilter.HasAnnotation ffann = (FeatureFilter.HasAnnotation)ff;
                ComparableTerm noteTerm;
                if (ffann.getKey() instanceof ComparableTerm) noteTerm = (ComparableTerm)ffann.getKey();
                else if (ffann.getKey() instanceof String) noteTerm =  RichObjectFactory.getDefaultOntology().getOrCreateTerm((String)ffann.getKey());
                else return null;
                return new BioSQLFeatureFilter.ByNoteTermOnly(noteTerm);
            }
            // ByAnnotation - attempt to convert the key to a ComparableTerm, and the value to a string, then wrap the
            // whole thing in a Note with rank 0 and try using ByNoteWithValue to retrieve.
            else if (ff instanceof FeatureFilter.ByAnnotation) {
                FeatureFilter.ByAnnotation ffann = (FeatureFilter.ByAnnotation)ff;
                if (!(ffann.getValue() instanceof String)) return null;
                String noteValue = (String)ffann.getValue();
                ComparableTerm noteTerm;
                Object key = ffann.getKey();
                if (key instanceof ComparableTerm) noteTerm = (ComparableTerm)key;
                else if (key instanceof String) noteTerm =  RichObjectFactory.getDefaultOntology().getOrCreateTerm((String)key);
                else return null;
                return new BioSQLFeatureFilter.ByNote(new SimpleNote(noteTerm,noteValue,0));
            }
            // OverlapsLocation - simple pass-through to OverlapsRichLocation.
            else if (ff instanceof FeatureFilter.OverlapsLocation) {
                FeatureFilter.OverlapsLocation ffloc = (FeatureFilter.OverlapsLocation)ff;
                return new BioSQLFeatureFilter.OverlapsRichLocation(RichLocation.Tools.enrich(ffloc.getLocation()));
            }
            // ShadowContainedByLocation - simple pass-through to ContainedByRichLocation, as we have no concept
            // of shadows within BioSQL so they are effectively the same thing.
            else if (ff instanceof FeatureFilter.ShadowContainedByLocation) {
                FeatureFilter.ShadowContainedByLocation ffloc = (FeatureFilter.ShadowContainedByLocation)ff;
                return new BioSQLFeatureFilter.ContainedByRichLocation(RichLocation.Tools.enrich(ffloc.getLocation()));
            }
            // Anything else we don't recognise? Return null!
            else {
                return null;
            }
        }
    }
    
    // Now for some useful filters.
    
    /**
     * All features are selected by this filter.
     */
    static final public BioSQLFeatureFilter all = new BioSQLAcceptAllFilter();
    
    /**
     * No features are selected by this filter.
     */
    static final public BioSQLFeatureFilter none = new BioSQLAcceptNoneFilter();
    
    /**
     * A filter for Hibernate-BioSQL filters to extend.
     */
    public abstract static class HibernateFeatureFilter implements BioSQLFeatureFilter {
        protected Method not;
        protected Method and;
        protected Method or;
        protected Method eq;
        protected Method le;
        protected Method ge;
        protected Method conjunction;
        protected Method disjunction;
        protected Method conjunctAdd;
        protected Method disjunctAdd;
        
        public HibernateFeatureFilter() {
            try {
                // Lazy load the Restrictions class from Hibernate.
                Class restrictions = Class.forName("org.hibernate.criterion.Restrictions");
                // Lazy load the Criterion class from Hibernate.
                Class criterion = Class.forName("org.hibernate.criterion.Criterion");
                // Lookup the methods
                this.not = restrictions.getMethod("not", new Class[]{criterion});
                this.and = restrictions.getMethod("and", new Class[]{criterion,criterion});
                this.or = restrictions.getMethod("or", new Class[]{criterion,criterion});
                this.eq = restrictions.getMethod("eq", new Class[]{String.class,Object.class});
                this.le = restrictions.getMethod("le", new Class[]{String.class,Object.class});
                this.ge = restrictions.getMethod("ge", new Class[]{String.class,Object.class});
                this.conjunction = restrictions.getMethod("conjunction", new Class[]{});
                this.disjunction = restrictions.getMethod("disjunction", new Class[]{});
                // Lazy load the Conjunction(Or)+Disjunction(And) class from Hibernate.
                Class conjunctClass = Class.forName("org.hibernate.criterion.Conjunction");
                Class disjunctClass = Class.forName("org.hibernate.criterion.Disjunction");
                // Lookup the methods
                this.conjunctAdd = conjunctClass.getMethod("add", new Class[]{criterion});
                this.disjunctAdd = disjunctClass.getMethod("add", new Class[]{criterion});
            } catch (ClassNotFoundException e) {
                throw new RuntimeException(e);
            } catch (NoSuchMethodException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            return Collections.EMPTY_MAP;
        }
    }
    
    /**
     *  A filter that returns all features not accepted by a child filter.
     *
     * @author Thomas Down
     * @author Matthew Pocock
     * @author Richard Holland
     * @since 1.5
     */
    public final static class Not extends HibernateFeatureFilter {
        static { WalkerFactory.getInstance().addTypeWithParent(Not.class); }
        
        BioSQLFeatureFilter child;
        
        public BioSQLFeatureFilter getChild() {
            return child;
        }
        
        public Not(BioSQLFeatureFilter child) {
            super();
            if (!(child instanceof BioSQLFeatureFilter))
                throw new BioRuntimeException("Cannot use non-BioSQLFeatureFilter instances with this class");
            this.child = child;
        }
        
        public boolean accept(Feature f) {
            return !(child.accept(f));
        }
        
        public Object asCriterion() {
            try {
                return this.not.invoke(null,new Object[]{child.asCriterion()});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            return child.criterionAliasMap();
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof Not) &&
                    (((Not) o).getChild().equals(this.getChild()));
        }
        
        public int hashCode() {
            return getChild().hashCode();
        }
        
        public String toString() {
            return "Not(" + child + ")";
        }
    }
    
    
    /**
     *  A filter that returns all features accepted by both child filter.
     *
     * @author Thomas Down
     * @author Matthew Pocock
     * @author Richard Holland
     * @since 1.5
     */
    public final static class And extends HibernateFeatureFilter {
        static { WalkerFactory.getInstance().addTypeWithParent(And.class); }
        
        BioSQLFeatureFilter c1, c2;
        
        public BioSQLFeatureFilter getChild1() {
            return c1;
        }
        
        public BioSQLFeatureFilter getChild2() {
            return c2;
        }
        
        public And(BioSQLFeatureFilter c1, BioSQLFeatureFilter c2) {
            super();
            if (!(c1 instanceof BioSQLFeatureFilter) || !(c2 instanceof BioSQLFeatureFilter))
                throw new BioRuntimeException("Cannot use non-BioSQLFeatureFilter instances with this class");
            this.c1 = c1;
            this.c2 = c2;
        }
        
        public boolean accept(Feature f) {
            return (c1.accept(f) && c2.accept(f));
        }
        
        public Object asCriterion() {
            try {
                return this.and.invoke(null,new Object[]{c1.asCriterion(),c2.asCriterion()});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            Map results = new HashMap();
            results.putAll(c1.criterionAliasMap());
            results.putAll(c2.criterionAliasMap());
            return results;
        }
        
        public boolean equals(Object o) {
            if(o instanceof BioSQLFeatureFilter) {
                return FilterUtils.areEqual(this, (FeatureFilter) o);
            } else {
                return false;
            }
        }
        
        public int hashCode() {
            return getChild1().hashCode() ^ getChild2().hashCode();
        }
        
        public String toString() {
            return "And(" + c1 + " , " + c2 + ")";
        }
    }
    
    /**
     *  A filter that returns all features accepted by at least one child filter.
     *
     * @author Thomas Down
     * @author Matthew Pocock
     * @author Richard Holland
     * @since 1.5
     */
    public final static class Or extends HibernateFeatureFilter {
        static { WalkerFactory.getInstance().addTypeWithParent(Or.class); }
        
        BioSQLFeatureFilter c1, c2;
        
        public BioSQLFeatureFilter getChild1() {
            return c1;
        }
        
        public BioSQLFeatureFilter getChild2() {
            return c2;
        }
        
        public Or(BioSQLFeatureFilter c1, BioSQLFeatureFilter c2) {
            super();
            if (!(c1 instanceof BioSQLFeatureFilter) || !(c2 instanceof BioSQLFeatureFilter))
                throw new BioRuntimeException("Cannot use non-BioSQLFeatureFilter instances with this class");
            this.c1 = c1;
            this.c2 = c2;
        }
        
        public boolean accept(Feature f) {
            return (c1.accept(f) || c2.accept(f));
        }
        
        public Object asCriterion() {
            try {
                return this.or.invoke(null,new Object[]{c1.asCriterion(),c2.asCriterion()});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            Map results = new HashMap();
            results.putAll(c1.criterionAliasMap());
            results.putAll(c2.criterionAliasMap());
            return results;
        }
        
        public boolean equals(Object o) {
            if(o instanceof BioSQLFeatureFilter) {
                return FilterUtils.areEqual(this, (FeatureFilter) o);
            } else {
                return false;
            }
        }
        
        public int hashCode() {
            return getChild1().hashCode() ^ getChild2().hashCode();
        }
        
        public String toString() {
            return "Or(" + c1 + " , " + c2 + ")";
        }
    }
    
    /**
     * Construct one of these to filter features by display name.
     *
     * @author Richard Holland
     * @since 1.5
     */
    final public static class ByName extends HibernateFeatureFilter {
        private String name;
        
        public String getName() {
            return name;
        }
        
        /**
         * Create a ByType filter that filters in all features with type fields
         * equal to type.
         *
         * @param name  the String to match type fields against
         */
        public ByName(String name) {
            super();
            if (name == null) {
                throw new NullPointerException("Name may not be null");
            }
            this.name = name;
        }
        
        /**
         * Returns true if the feature has a matching type property.
         */
        public boolean accept(Feature f) {
            if (f instanceof RichFeature) {
                return name.equals(((RichFeature)f).getName());
            }
            return false;
        }
        
        public Object asCriterion() {
            try {
                return this.eq.invoke(null,new Object[]{"name",name});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof ByName) &&
                    (((ByName) o).getName().equals(this.getName()));
        }
        
        public int hashCode() {
            return getName().hashCode();
        }
        
        public String toString() {
            return "ByName(" + name + ")";
        }
    }
    
    /**
     * Construct one of these to filter features by rank.
     *
     * @author Richard Holland
     * @since 1.5
     */
    final public static class ByRank extends HibernateFeatureFilter {
        private int rank;
        
        public int getRank() {
            return rank;
        }
        
        /**
         * Create a Rank filter that filters in all features with rank fields
         * equal to rank.
         *
         * @param rank  the rank to match type fields against
         */
        public ByRank(int rank) {
            super();
            this.rank = rank;
        }
        
        /**
         * Returns true if the feature has a matching type property.
         */
        public boolean accept(Feature f) {
            if (f instanceof RichFeature) {
                return rank==((RichFeature)f).getRank();
            }
            return false;
        }
        
        public Object asCriterion() {
            try {
                return this.eq.invoke(null,new Object[]{"rank",new Integer(rank)});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof ByRank) &&
                    (((ByRank) o).getRank() == this.getRank());
        }
        
        public int hashCode() {
            return rank;
        }
        
        public String toString() {
            return "ByRank(" + rank + ")";
        }
    }
    
    /**
     * Construct one of these to filter features by type.
     *
     * @author Matthew Pocock
     * @author Richard Holland
     * @since 1.5
     */
    final public static class ByTypeTerm extends HibernateFeatureFilter {
        private ComparableTerm typeTerm;
        
        public ComparableTerm getTypeTerm() {
            return typeTerm;
        }
        
        /**
         * Create a ByTypeTerm filter that filters in all features with typeTerm fields
         * equal to typeTerm.
         *
         * @param typeTerm  the Term to match typeTerm fields against
         */
        public ByTypeTerm(ComparableTerm typeTerm) {
            super();
            if (typeTerm == null) {
                throw new NullPointerException("Type may not be null");
            }
            this.typeTerm = typeTerm;
        }
        
        /**
         * Returns true if the feature has a matching type property.
         */
        public boolean accept(Feature f) {
            return typeTerm.equals(f.getTypeTerm());
        }
        
        public Object asCriterion() {
            try {
                return this.eq.invoke(null,new Object[]{"typeTerm",typeTerm});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof ByTypeTerm) &&
                    (((ByTypeTerm) o).getTypeTerm().equals(this.getTypeTerm()));
        }
        
        public int hashCode() {
            return getTypeTerm().hashCode();
        }
        
        public String toString() {
            return "ByTypeTerm(" + typeTerm + ")";
        }
    }
    
    
    /**
     * Construct one of these to filter features by source.
     *
     * @author Matthew Pocock
     * @author Richard Holland
     * @since 1.5
     */
    final public static class BySourceTerm extends HibernateFeatureFilter {
        private ComparableTerm sourceTerm;
        
        public ComparableTerm getSourceTerm() {
            return sourceTerm;
        }
        
        /**
         * Create a BySourceTerm filter that filters in all features with sourceTerm fields
         * equal to source.
         *
         * @param sourceTerm  the Term to match sourceTerm fields against
         */
        public BySourceTerm(ComparableTerm sourceTerm) {
            super();
            if (sourceTerm == null) {
                throw new NullPointerException("Source may not be null");
            }
            this.sourceTerm = sourceTerm;
        }
        
        /**
         * Returns true if the feature has a matching source property.
         */
        public boolean accept(Feature f) {
            return sourceTerm.equals(f.getSourceTerm());
        }
        
        public Object asCriterion() {
            try {
                return this.eq.invoke(null,new Object[]{"sourceTerm",sourceTerm});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof BySourceTerm) &&
                    (((BySourceTerm) o).getSourceTerm().equals(this.getSourceTerm()));
        }
        
        public int hashCode() {
            return getSourceTerm().hashCode();
        }
        
        public String toString() {
            return "BySourceTerm(" + sourceTerm + ")";
        }
    }
    
    /**
     * Construct one of these to filter features by type (name only - parent ontology
     * is ignored).
     *
     * @author Richard Holland
     * @since 1.5
     */
    final public static class ByTypeTermName extends HibernateFeatureFilter {
        private String typeTermName;
        
        public String getTypeTermName() {
            return typeTermName;
        }
        
        /**
         * Create a ByTypeTermName filter that filters in all features with typeTerm fields
         * having name equal to typeTermName.
         *
         * @param typeTermName  the Term to match typeTermName fields against
         */
        public ByTypeTermName(String typeTermName) {
            super();
            if (typeTermName == null) {
                throw new NullPointerException("Type name may not be null");
            }
            this.typeTermName = typeTermName;
        }
        
        /**
         * Returns true if the feature has a matching type property.
         */
        public boolean accept(Feature f) {
            return typeTermName.equals(f.getTypeTerm().getName());
        }
        
        public Object asCriterion() {
            try {
                return this.eq.invoke(null,new Object[]{"tt.name",typeTermName});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            Map results = new HashMap();
            results.put("typeTerm","tt");
            return results;
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof ByTypeTermName) &&
                    (((ByTypeTermName) o).getTypeTermName().equals(this.getTypeTermName()));
        }
        
        public int hashCode() {
            return getTypeTermName().hashCode();
        }
        
        public String toString() {
            return "ByTypeTermName(" + typeTermName + ")";
        }
    }
    
    
    /**
     * Construct one of these to filter features by source (name only - parent ontology
     * is ignored).
     *
     * @author Richard Holland
     * @since 1.5
     */
    final public static class BySourceTermName extends HibernateFeatureFilter {
        private String sourceTermName;
        
        public String getSourceTermName() {
            return sourceTermName;
        }
        
        /**
         * Create a BySourceTerm filter that filters in all features with sourceTerm fields
         * having name equal to sourceTermName.
         *
         * @param sourceTermName  the name of the Term to match sourceTerm fields against
         */
        public BySourceTermName(String sourceTermName) {
            super();
            if (sourceTermName == null) {
                throw new NullPointerException("Source name may not be null");
            }
            this.sourceTermName = sourceTermName;
        }
        
        /**
         * Returns true if the feature has a matching source property.
         */
        public boolean accept(Feature f) {
            return sourceTermName.equals(f.getSourceTerm().getName());
        }
        
        public Object asCriterion() {
            try {
                return this.eq.invoke(null,new Object[]{"st.name",sourceTermName});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            Map results = new HashMap();
            results.put("sourceTerm","st");
            return results;
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof BySourceTermName) &&
                    (((BySourceTermName) o).getSourceTermName().equals(this.getSourceTermName()));
        }
        
        public int hashCode() {
            return getSourceTermName().hashCode();
        }
        
        public String toString() {
            return "BySourceTermName(" + sourceTermName + ")";
        }
    }
    
    /**
     * Accept features that reside on a sequence with a particular name.
     *
     * @author Matthew Pocock
     * @author Richard Holland
     * @since 1.5
     */
    public final static class BySequenceName extends HibernateFeatureFilter {
        private String seqName;
        
        public BySequenceName(String seqName) {
            super();
            this.seqName = seqName;
        }
        
        public String getSequenceName() {
            return seqName;
        }
        
        public boolean accept(Feature f) {
            return f.getSequence().getName().equals(seqName);
        }
        
        public Object asCriterion() {
            try {
                return this.eq.invoke(null,new Object[]{"p.name",seqName});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            Map results = new HashMap();
            results.put("parent","p");
            return results;
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof BySequenceName) &&
                    ((BySequenceName) o).getSequenceName().equals(seqName);
        }
        
        public int hashCode() {
            return seqName.hashCode();
        }
    }
    
    /**
     * A filter that returns all features contained within a location. Contained means
     * that a feature is entirely within, on the same strand and on the same sequence
     * as any single member of the flattened query location.
     *
     * @author Matthew Pocock
     * @author Richard Holland
     * @since 1.5
     */
    public final static class ContainedByRichLocation extends HibernateFeatureFilter {
        private RichLocation loc;
        
        public RichLocation getRichLocation() {
            return loc;
        }
        
        /**
         * Creates a filter that returns everything contained within loc.
         *
         * @param loc  the location that will contain the accepted features
         */
        public ContainedByRichLocation(RichLocation loc) {
            super();
            if (loc == null) {
                throw new NullPointerException("Loc may not be null");
            }
            this.loc = loc;
        }
        
        /**
         * Returns true if the feature is within this filter's location.
         */
        public boolean accept(Feature f) {
            return loc.contains(f.getLocation());
        }
        
        public Object asCriterion() {
            try {
                // Conjunction of criteria for each member of the query location.
                Collection members = RichLocation.Tools.flatten(loc);
                // some combo of Tools.flatten(loc), min(loc.start,feat.start) and min(loc.end,feat.end)
                Object parentConjunct = this.conjunction.invoke((Object[])null,(Object[])null);
                for (Iterator i = members.iterator(); i.hasNext(); ) {
                    RichLocation loc = (RichLocation)i.next();
                    Object childDisjunct = this.disjunction.invoke((Object[])null,(Object[])null);
                    // for each member, find features that have start>=member.start,
                    // end<=member.end and strand=member.strand and crossref=member.crossref
                    this.disjunctAdd.invoke(childDisjunct,new Object[]{this.eq.invoke(null, new Object[]{"l.strandNum",new Integer(loc.getStrand().intValue())})});
                    this.disjunctAdd.invoke(childDisjunct,new Object[]{this.eq.invoke(null, new Object[]{"l.crossRef",loc.getCrossRef()})});
                    this.disjunctAdd.invoke(childDisjunct,new Object[]{this.ge.invoke(null, new Object[]{"l.min",new Integer(loc.getMin())})});
                    this.disjunctAdd.invoke(childDisjunct,new Object[]{this.le.invoke(null, new Object[]{"l.max",new Integer(loc.getMax())})});
                    // add the member to the set of restrictions
                    this.conjunctAdd.invoke(parentConjunct,new Object[]{childDisjunct});
                }
                return parentConjunct;
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            Map results = new HashMap();
            results.put("locationSet","l");
            return results;
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof ContainedByRichLocation) &&
                    (((ContainedByRichLocation) o).getRichLocation().equals(this.getRichLocation()));
        }
        
        public int hashCode() {
            return getRichLocation().hashCode();
        }
        
        public String toString() {
            return "ContainedBy(" + loc + ")";
        }
    }
    
    /**
     * A filter that returns all features having locations on a given strand. They
     * may actually have features on other strands too, of course.
     *
     * @author Richard Holland
     * @since 1.5
     */
    public final static class ByStrand extends HibernateFeatureFilter {
        private Strand str;
        
        public Strand getStrand() {
            return str;
        }
        
        /**
         * Creates a filter that returns everything on strand str.
         *
         * @param str  the strand that will contain the accepted features
         */
        public ByStrand(Strand str) {
            super();
            if (str == null) {
                throw new NullPointerException("Strand may not be null");
            }
            this.str = str;
        }
        
        /**
         * Returns true if the feature overlaps this filter's location.
         */
        public boolean accept(Feature f) {
            if (f instanceof RichFeature) {
                RichFeature rf = (RichFeature)f;
                for (Iterator i = rf.getLocation().blockIterator(); i.hasNext(); ) {
                    RichLocation l = (RichLocation)i.next();
                    if (l.getStrand().equals(str)) return true;
                }
            }
            return false;
        }
        
        public Object asCriterion() {
            try {
                // any location on the feature with a matching strand?
                return this.eq.invoke(null,new Object[]{"l.strandNum",new Integer(str.intValue())});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            Map results = new HashMap();
            results.put("locationSet","l");
            return results;
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof ByStrand) &&
                    (((ByStrand) o).getStrand().equals(this.getStrand()));
        }
        
        public int hashCode() {
            return getStrand().hashCode();
        }
        
        public String toString() {
            return "ByStrand(" + str + ")";
        }
    }
    
    /**
     * A filter that returns all features overlapping a location. Overlaps means
     * that a feature includes part of, on the same strand and on the same sequence
     * any single member of the flattened query location.
     *
     * @author Matthew Pocock
     * @author Richard Holland
     * @since 1.5
     */
    public final static class OverlapsRichLocation extends HibernateFeatureFilter {
        private RichLocation loc;
        
        public RichLocation getRichLocation() {
            return loc;
        }
        
        /**
         * Creates a filter that returns everything overlapping loc.
         *
         * @param loc  the location that will overlap the accepted features
         */
        public OverlapsRichLocation(RichLocation loc) {
            super();
            if (loc == null) {
                throw new NullPointerException("Loc may not be null");
            }
            this.loc = loc;
        }
        
        /**
         * Returns true if the feature overlaps this filter's location.
         */
        public boolean accept(Feature f) {
            return loc.overlaps(f.getLocation());
        }
        
        public Object asCriterion() {
            try {
                // Conjunction of criteria for each member of the query location.
                Collection members = RichLocation.Tools.flatten(loc);
                // some combo of Tools.flatten(loc), min(loc.start,feat.start) and min(loc.end,feat.end)
                Object parentConjunct = this.conjunction.invoke((Object[])null,(Object[])null);
                for (Iterator i = members.iterator(); i.hasNext(); ) {
                    RichLocation loc = (RichLocation)i.next();
                    Object childDisjunct = this.disjunction.invoke((Object[])null,(Object[])null);
                    // for each member, find features that have start<=member.end,  end>=member.start,
                    // strand=member.strand and crossref=member.crossref
                    this.disjunctAdd.invoke(childDisjunct,new Object[]{this.eq.invoke(null, new Object[]{"l.strandNum",new Integer(loc.getStrand().intValue())})});
                    this.disjunctAdd.invoke(childDisjunct,new Object[]{this.eq.invoke(null, new Object[]{"l.crossRef",loc.getCrossRef()})});
                    this.disjunctAdd.invoke(childDisjunct,new Object[]{this.ge.invoke(null, new Object[]{"l.max",new Integer(loc.getMin())})});
                    this.disjunctAdd.invoke(childDisjunct,new Object[]{this.le.invoke(null, new Object[]{"l.min",new Integer(loc.getMax())})});
                    // add the member to the set of restrictions
                    this.conjunctAdd.invoke(parentConjunct,new Object[]{childDisjunct});
                }
                return parentConjunct;
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            Map results = new HashMap();
            results.put("locationSet","l");
            return results;
        }
        
        public boolean equals(Object o) {
            return
                    (o instanceof OverlapsRichLocation) &&
                    (((OverlapsRichLocation) o).getRichLocation().equals(this.getRichLocation()));
        }
        
        public int hashCode() {
            return getRichLocation().hashCode();
        }
        
        public String toString() {
            return "Overlaps(" + loc + ")";
        }
    }
    
    /**
     * A filter that returns all features that have the given note, and
     * the value and rank is checked as well.
     *
     * @author Richard Holland
     * @author George Waldon
     * @since 1.5
     */
    public final static class ByNote extends HibernateFeatureFilter {
        private Note note;
        
        public ByNote(Note note) {
            super();
            this.note = note;
        }
        
        public Note getNote() {
            return note;
        }
        
        public boolean accept(Feature f) {
            if (f instanceof RichFeature) {
                RichAnnotation ra = ((RichFeature)f).getRichAnnotation();
                try {
                    Note n = ra.getNote(note);
                    return (n.getValue()==note.getValue()) || (n.getValue()!=null && note.getValue()!=null && n.getValue().equals(note.getValue()));
                } catch (NoSuchElementException e) {
                    return false;
                }
            }
            return false;
        }
        
        public Object asCriterion() {
            try {
                Object conjunct = this.conjunction.invoke((Object[])null,(Object[])null);
                this.conjunctAdd.invoke(conjunct,new Object[]{this.eq.invoke(null, new Object[]{"n.term",note.getTerm()})});
                this.conjunctAdd.invoke(conjunct,new Object[]{this.eq.invoke(null, new Object[]{"n.value",note.getValue()})});
                this.conjunctAdd.invoke(conjunct,new Object[]{this.eq.invoke(null, new Object[]{"n.rank",new Integer(note.getRank())})});
                return conjunct;
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            Map results = new HashMap();
            results.put("noteSet","n");
            return results;
        }
        
        public boolean equals(Object o) {
            if(o instanceof ByNote) {
                ByNote that = (ByNote) o;
                return this.getNote() == that.getNote();
            }
            
            return false;
        }
        
        public int hashCode() {
            return getNote().hashCode();
        }
        
        public String toString() {
            return "ByNote {" + note + "}";
        }
    }
    
    /**
     * A filter that returns all features that have a note with the given term. The value
     * and rank is not checked.
     *
     * @author Richard Holland
     * @author George Waldon
     * @since 1.5
     */
    public final static class ByNoteTermOnly extends HibernateFeatureFilter {
        private ComparableTerm term;
        
        public ByNoteTermOnly(ComparableTerm term) {
            super();
            this.term = term;
        }
        
        public ComparableTerm getTerm() {
            return term;
        }
        
        public boolean accept(Feature f) {
            if (f instanceof RichFeature) {
                RichAnnotation ra = ((RichFeature)f).getRichAnnotation();
                try {
                    for (Iterator i = ra.getNoteSet().iterator(); i.hasNext(); ) if (((Note)i.next()).getTerm().equals(term)) return true;
                } catch (NoSuchElementException e) {
                    return false;
                }
            }
            return false;
        }
        
        public Object asCriterion() {
            try {
                return this.eq.invoke(null, new Object[]{"n.term",term});
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        
        public Map criterionAliasMap() {
            Map results = new HashMap();
            results.put("noteSet","n");
            return results;
        }
        
        public boolean equals(Object o) {
            if(o instanceof ByNoteTermOnly) {
                ByNoteTermOnly that = (ByNoteTermOnly) o;
                return this.getTerm() == that.getTerm();
            }
            
            return false;
        }
        
        public int hashCode() {
            return getTerm().hashCode();
        }
        
        public String toString() {
            return "ByNoteTermOnly {" + term + "}";
        }
    }
}

