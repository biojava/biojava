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

package org.biojava.bio.seq;

import junit.framework.TestCase;

import org.biojava.bio.AnnotationType;
import org.biojava.bio.CardinalityConstraint;
import org.biojava.bio.PropertyConstraint;
import org.biojava.bio.seq.homol.HomologyFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.RangeLocation;

/**
 * Tests for FeatureFilters.  Currently concentrating on
 * properSubset and disjunction.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.2
 */
public class FilterUtilsTest extends TestCase
{
    protected FeatureFilter all_and_all;
    protected FeatureFilter all_and_none;
    protected FeatureFilter none_and_all;
    protected FeatureFilter none_and_none;
    
    protected FeatureFilter all_or_all;
    protected FeatureFilter all_or_none;
    protected FeatureFilter none_or_all;
    protected FeatureFilter none_or_none;
    
    protected FeatureFilter tf1;
    protected FeatureFilter tf2;
    protected FeatureFilter tf3;

    protected FeatureFilter pf1;
    protected FeatureFilter pf2;
    protected FeatureFilter pf3;
    protected FeatureFilter pf4;
    protected FeatureFilter pf5;
    protected FeatureFilter pf6;
   
    protected FeatureFilter cf_StrandedFeature;
    protected FeatureFilter cf_ComponentFeature;
    protected FeatureFilter cf_HomologyFeature;

    protected FeatureFilter parent_cf_ComponentFeature;
    protected FeatureFilter ancestor_cf_ComponentFeature;
    protected FeatureFilter not_parent_cf_ComponentFeature;
    protected FeatureFilter not_ancestor_cf_ComponentFeature;

    protected FeatureFilter ntf1;
    protected FeatureFilter ntf2;

    protected FeatureFilter olf1;
    protected FeatureFilter olf2;
    protected FeatureFilter olf3;
    protected FeatureFilter olf4;
    
    protected FeatureFilter clf1;
    protected FeatureFilter clf2;
    protected FeatureFilter clf3;
    protected FeatureFilter clf4;

    protected FeatureFilter tf1_or_tf2;
    protected FeatureFilter tf2_or_tf3;
    protected FeatureFilter tf1_or_tf3;
    protected FeatureFilter tf1_or_tf2_or_tf3;

    protected FeatureFilter pf1_and_pf2;
    protected FeatureFilter pf2_and_pf3;
    protected FeatureFilter pf1_and_pf3;
    protected FeatureFilter pf1_and_pf2_and_pf3;

    protected FeatureFilter pf1_and_tf1;
    protected FeatureFilter pf1_and_pf2_and_tf1;

    protected FeatureFilter pf1_AND_tf1_or_tf2;
    protected FeatureFilter pf1_and_pf2_OR_tf1;

    public FilterUtilsTest(String name) {
	super(name);
    }

    protected void setUp() throws Exception {
      //
      // pure and logicals
      //
      all_and_all = new FeatureFilter.And(FeatureFilter.all, FeatureFilter.all);
      all_and_none = new FeatureFilter.And(FeatureFilter.all, FeatureFilter.none);
      none_and_all = new FeatureFilter.And(FeatureFilter.none, FeatureFilter.all);
      none_and_none = new FeatureFilter.And(FeatureFilter.none, FeatureFilter.none);
      
      //
      // pure or logicals
      //
      all_or_all = new FeatureFilter.Or(FeatureFilter.all, FeatureFilter.all);
      all_or_none = new FeatureFilter.Or(FeatureFilter.all, FeatureFilter.none);
      none_or_all = new FeatureFilter.Or(FeatureFilter.none, FeatureFilter.all);
      none_or_none = new FeatureFilter.Or(FeatureFilter.none, FeatureFilter.none);

  //
	// Three type filters (opaque but mutually disjoint).
	//

	tf1 = new FeatureFilter.ByType("hello");
	tf2 = new FeatureFilter.ByType("goodbye");
	tf3 = new FeatureFilter.ByType("moo");

	//
	// Three annotation-property filters (opaque, non-disjoint)
	//

	pf1 = new FeatureFilter.HasAnnotation("foo");
	pf2 = new FeatureFilter.HasAnnotation("bar");
	pf3 = new FeatureFilter.HasAnnotation("baz");

  //
  // Three annotation-property filters with values
  //
  
  pf4 = new FeatureFilter.ByAnnotation("foo", "fish");
  pf5 = new FeatureFilter.ByAnnotation("foo", "cat");
  pf6 = new FeatureFilter.ByAnnotation("bar", "fish");
  
	//
	// Class filters
	//
	
	cf_StrandedFeature = new FeatureFilter.ByClass(StrandedFeature.class);
	cf_ComponentFeature = new FeatureFilter.ByClass(ComponentFeature.class);
	cf_HomologyFeature = new FeatureFilter.ByClass(HomologyFeature.class);

	//
	// NOTed filters.
	//

	ntf1 = new FeatureFilter.Not(tf1);
	ntf2 = new FeatureFilter.Not(tf2);

	//
	// Various locations
	//

	Location l1 = new RangeLocation(1000, 2000);
	Location l2 = new RangeLocation(10000, 11000);
	Location l3 = new PointLocation(10500);
	Location l4 = new RangeLocation(10700, 11700);

	olf1 = new FeatureFilter.OverlapsLocation(l1);
	olf2 = new FeatureFilter.OverlapsLocation(l2);
	olf3 = new FeatureFilter.OverlapsLocation(l3);
	olf4 = new FeatureFilter.OverlapsLocation(l4);

	clf1 = new FeatureFilter.ContainedByLocation(l1);
	clf2 = new FeatureFilter.ContainedByLocation(l2);
	clf3 = new FeatureFilter.ContainedByLocation(l3);
	clf4 = new FeatureFilter.ContainedByLocation(l4);

	//
	// ANDs and ORs
	//

	tf1_or_tf2 = new FeatureFilter.Or(tf1, tf2);
	tf2_or_tf3 = new FeatureFilter.Or(tf2, tf3);
	tf1_or_tf3 = new FeatureFilter.Or(tf1, tf3);
	tf1_or_tf2_or_tf3 = new FeatureFilter.Or(tf1_or_tf2, tf3);

	pf1_and_pf2 = new FeatureFilter.And(pf1, pf2);
	pf2_and_pf3 = new FeatureFilter.And(pf2, pf3);
	pf1_and_pf3 = new FeatureFilter.And(pf1, pf3);
	pf1_and_pf2_and_pf3 = new FeatureFilter.And(pf1_and_pf2, pf3);
	pf1_and_tf1 = new FeatureFilter.And(pf1, tf1);
	pf1_and_pf2_and_tf1 = new FeatureFilter.And(pf1_and_pf2, tf1);

	pf1_AND_tf1_or_tf2 = new FeatureFilter.And(pf1, tf1_or_tf2);
	pf1_and_pf2_OR_tf1 = new FeatureFilter.Or(pf1_and_pf2, tf1);

	//
	// Parent and Ancestor filters.
	//
	
	parent_cf_ComponentFeature = new FeatureFilter.ByParent(cf_ComponentFeature);
	ancestor_cf_ComponentFeature = new FeatureFilter.ByAncestor(cf_ComponentFeature);
	not_parent_cf_ComponentFeature = new FeatureFilter.Not(parent_cf_ComponentFeature);
	not_ancestor_cf_ComponentFeature = new FeatureFilter.Not(ancestor_cf_ComponentFeature);
	
    }

    public void testTypes() throws Exception {
	assertTrue(! FilterUtils.areProperSubset(tf1, tf2));
	assertTrue(FilterUtils.areProperSubset(tf1, tf1));

	assertTrue(FilterUtils.areDisjoint(tf1, tf2));
	assertTrue(! FilterUtils.areDisjoint(tf1, tf1));
    }

    public void testByClass() throws Exception {
	assertTrue(FilterUtils.areProperSubset(cf_HomologyFeature, cf_StrandedFeature));
	assertTrue(FilterUtils.areProperSubset(cf_ComponentFeature, cf_StrandedFeature));
	assertTrue(FilterUtils.areDisjoint(cf_ComponentFeature, cf_HomologyFeature));
    }

    public void testLocation() throws Exception {
	// Simple comparison of locations.

	assertTrue(FilterUtils.areProperSubset(olf1, olf1));
	assertTrue(! FilterUtils.areProperSubset(olf1, olf2));
	assertTrue(! FilterUtils.areDisjoint(olf1, olf1));
	
	// Assymetry between containment and overlapping

	assertTrue(FilterUtils.areProperSubset(clf1, olf1));
	assertTrue(! FilterUtils.areProperSubset(olf1, clf1));

	assertTrue(FilterUtils.areDisjoint(olf1, clf2));
	assertTrue(FilterUtils.areDisjoint(clf1, olf2));
    }

    public void testNot() throws Exception {
	assertTrue(FilterUtils.areProperSubset(ntf1, ntf1));
	assertTrue(! FilterUtils.areDisjoint(ntf1, ntf1));

	assertTrue(FilterUtils.areDisjoint(tf1, ntf1));
	assertTrue(! FilterUtils.areProperSubset(tf1, ntf1));

	assertTrue(! FilterUtils.areDisjoint(tf2, ntf1));
	assertTrue(FilterUtils.areProperSubset(tf2, ntf1));
    }

    public void testOr() throws Exception {
	assertTrue(FilterUtils.areProperSubset(tf1, tf1_or_tf2));
	assertTrue(! FilterUtils.areDisjoint(tf1, tf1_or_tf2));
	assertTrue(FilterUtils.areDisjoint(tf3, tf1_or_tf2));
	
	assertTrue(! FilterUtils.areProperSubset(tf1_or_tf2,
						 tf2_or_tf3));
	assertTrue(! FilterUtils.areDisjoint(tf1_or_tf2,
					     tf2_or_tf3));

	assertTrue(FilterUtils.areProperSubset(tf1_or_tf2,
					       tf1_or_tf2_or_tf3));
	assertTrue(FilterUtils.areProperSubset(tf2_or_tf3,
					       tf1_or_tf2_or_tf3));
	assertTrue(FilterUtils.areProperSubset(tf1_or_tf3,
					       tf1_or_tf2_or_tf3));
    }

    public void testAnd() throws Exception {
      //all_and_all vs all
      assertTrue("are subset: " + all_and_all + ", " + FeatureFilter.all, FilterUtils.areProperSubset(all_and_all, FeatureFilter.all));
      assertTrue("not disjoint: " + all_and_all + ", " + FeatureFilter.all, !FilterUtils.areDisjoint(all_and_all, FeatureFilter.all));
      assertTrue("are subset: " + FeatureFilter.all + ", " + all_and_all, FilterUtils.areProperSubset(FeatureFilter.all, all_and_all));
      assertTrue("not disjoint: " + FeatureFilter.all + ", " + all_and_all, !FilterUtils.areDisjoint(FeatureFilter.all, all_and_all));

      // all_and_all vs none
      assertTrue("not subset: " + all_and_all + ", " + FeatureFilter.none, !FilterUtils.areProperSubset(all_and_all, FeatureFilter.none));
      assertTrue("are disjoint: " + all_and_all + ", " + FeatureFilter.none, FilterUtils.areDisjoint(all_and_all, FeatureFilter.none));
      assertTrue("are subset: " + FeatureFilter.none + ", " + all_and_all, FilterUtils.areProperSubset(FeatureFilter.none, all_and_all));
      assertTrue("are disjoint: " + FeatureFilter.none + ", " + all_and_all, FilterUtils.areDisjoint(FeatureFilter.none, all_and_all));
      
      // all_and_none vs all
      assertTrue("are subset: " + all_and_none + ", " + FeatureFilter.all, FilterUtils.areProperSubset(all_and_none, FeatureFilter.all));
      assertTrue("are disjoint: " + all_and_none + ", " + FeatureFilter.all, FilterUtils.areDisjoint(all_and_none, FeatureFilter.all));
      assertTrue("not subset: " + FeatureFilter.all + ", " + all_and_none, !FilterUtils.areProperSubset(FeatureFilter.all, all_and_none));
      assertTrue("are disjoint: " + FeatureFilter.all + ", " + all_and_none, FilterUtils.areDisjoint(FeatureFilter.all, all_and_none));
      
      // all_and_none vs none
      assertTrue("are subset: " + all_and_none + ", " + FeatureFilter.none, FilterUtils.areProperSubset(all_and_none, FeatureFilter.none));
      assertTrue("are disjoint: " + all_and_none + ", " + FeatureFilter.none, FilterUtils.areDisjoint(all_and_none, FeatureFilter.none));
      assertTrue("are subset: " + FeatureFilter.none + ", " + all_and_none, FilterUtils.areProperSubset(FeatureFilter.none, all_and_none));
      assertTrue("are disjoint: " + FeatureFilter.none + ", " + all_and_none, FilterUtils.areDisjoint(FeatureFilter.none, all_and_none));

      // none_and_all vs all
      assertTrue("are subset: " + none_and_all + ", " + FeatureFilter.all, FilterUtils.areProperSubset(none_and_all, FeatureFilter.all));
      assertTrue("are disjoint: " + none_and_all + ", " + FeatureFilter.all, FilterUtils.areDisjoint(none_and_all, FeatureFilter.all));
      assertTrue("not subset: " + FeatureFilter.all + ", " + none_and_all, !FilterUtils.areProperSubset(FeatureFilter.all, none_and_all));
      assertTrue("are disjoint: " + FeatureFilter.all + ", " + none_and_all, FilterUtils.areDisjoint(FeatureFilter.all, none_and_all));
      
      // none_and_all vs none
      assertTrue("are subset: " + none_and_all + ", " + FeatureFilter.none, FilterUtils.areProperSubset(none_and_all, FeatureFilter.none));
      assertTrue("are disjoint: " + none_and_all + ", " + FeatureFilter.none, FilterUtils.areDisjoint(none_and_all, FeatureFilter.none));
      assertTrue("are subset: " + FeatureFilter.none + ", " + none_and_all, FilterUtils.areProperSubset(FeatureFilter.none, none_and_all));
      assertTrue("are disjoint: " + FeatureFilter.none + ", " + none_and_all, FilterUtils.areDisjoint(FeatureFilter.none, none_and_all));
      
      // none_and_none vs all
      assertTrue("are subset: " + none_and_none + ", " + FeatureFilter.all, FilterUtils.areProperSubset(none_and_none, FeatureFilter.all));
      assertTrue("are disjoint: " + none_and_none + ", " + FeatureFilter.all, FilterUtils.areDisjoint(none_and_none, FeatureFilter.all));
      assertTrue("not subset: " + FeatureFilter.all + ", " + none_and_none, !FilterUtils.areProperSubset(FeatureFilter.all, none_and_none));
      assertTrue("are disjoint: " + FeatureFilter.all + ", " + none_and_none, FilterUtils.areDisjoint(FeatureFilter.all, none_and_none));
      
      // none_and_none vs none
      assertTrue("are subset: " + none_and_none + ", " + FeatureFilter.none, FilterUtils.areProperSubset(none_and_none, FeatureFilter.none));
      assertTrue("are disjoint: " + none_and_none + ", " + FeatureFilter.none, FilterUtils.areDisjoint(none_and_none, FeatureFilter.none));
      assertTrue("are subset: " + FeatureFilter.none + ", " + none_and_none, FilterUtils.areProperSubset(FeatureFilter.none, none_and_none));
      assertTrue("are disjoint: " + FeatureFilter.none + ", " + none_and_none, FilterUtils.areDisjoint(FeatureFilter.none, none_and_none));
      
      // pf1_and_pf2 vs pf1, pf2
      assertTrue("are subset: " + pf1_and_pf2 + ", " + pf1, FilterUtils.areProperSubset(pf1_and_pf2, pf1));
      assertTrue("not disjoint: " + pf1_and_pf2 + ", " + pf1, !FilterUtils.areDisjoint(pf1_and_pf2, pf1));
      assertTrue("are subset: " + pf1_and_pf2 + ", " + pf2, FilterUtils.areProperSubset(pf1_and_pf2, pf2));
      assertTrue("not disjoint: " + pf1_and_pf2 + ", " + pf2, !FilterUtils.areDisjoint(pf1_and_pf2, pf2));
      
      assertTrue("not subset: " + pf1_and_pf2 + ", " + pf3, !FilterUtils.areProperSubset(pf1_and_pf2, pf3));
      assertTrue("not subset: " + pf2_and_pf3 + ", " + pf1_and_pf2, !FilterUtils.areProperSubset(pf2_and_pf3, pf1_and_pf2));
      assertTrue("are subset: " + pf1_and_pf2_and_pf3 + ", " + pf1_and_pf2, FilterUtils.areProperSubset(pf1_and_pf2_and_pf3, pf1_and_pf2));
      assertTrue("are subset: " + pf1_and_pf2_and_pf3 + ", " + pf2_and_pf3, FilterUtils.areProperSubset(pf1_and_pf2_and_pf3, pf2_and_pf3));
      
      assertTrue("are subset: " + pf1_and_pf2_and_tf1 + ", " + pf1_and_pf2, FilterUtils.areProperSubset(pf1_and_pf2_and_tf1, pf1_and_pf2));
      assertTrue("are subset: " + pf1_and_pf2_and_tf1 + ", " + pf1_and_tf1, FilterUtils.areProperSubset(pf1_and_pf2_and_tf1, pf1_and_tf1));
      
      assertTrue("are disjoint: " + pf1_and_tf1 + ", " + tf2, FilterUtils.areDisjoint(pf1_and_tf1, tf2));
      assertTrue("are disjoint: " + pf1_and_pf2_and_tf1 + ", " + tf2, FilterUtils.areDisjoint(pf1_and_pf2_and_tf1, tf2));	   
    }

    public void testAndOr() throws Exception {
	assertTrue(FilterUtils.areProperSubset(pf1_and_tf1, pf1_AND_tf1_or_tf2));
	assertTrue(FilterUtils.areProperSubset(pf1_and_pf2_and_tf1, pf1_and_pf2_OR_tf1));
    }

    public void testParentAncestor() throws Exception {
	assertTrue(FilterUtils.areProperSubset(parent_cf_ComponentFeature, ancestor_cf_ComponentFeature));
	assertTrue(! FilterUtils.areProperSubset(ancestor_cf_ComponentFeature, parent_cf_ComponentFeature));
	assertTrue(FilterUtils.areDisjoint(ancestor_cf_ComponentFeature, not_ancestor_cf_ComponentFeature));
	assertTrue(FilterUtils.areDisjoint(not_ancestor_cf_ComponentFeature, parent_cf_ComponentFeature));
    }
    
    public void testHasProperty() throws Exception {
      assertTrue("not disjoint: " + pf1 + ", " + pf2, !FilterUtils.areDisjoint(pf1, pf2));
      assertTrue("not disjoint: " + pf2 + ", " + pf3, !FilterUtils.areDisjoint(pf2, pf3));
      assertTrue("not disjoint: " + pf1 + ", " + pf3, !FilterUtils.areDisjoint(pf1, pf3));
      
      assertTrue("not subset: " + pf1 + ", " + pf2, !FilterUtils.areProperSubset(pf1, pf2));
      assertTrue("not subset: " + pf1 + ", " + pf3, !FilterUtils.areProperSubset(pf1, pf3));
      assertTrue("not subset: " + pf2 + ", " + pf3, !FilterUtils.areProperSubset(pf2, pf3));

      assertTrue("are subset: " + pf1 + ", " + pf1, FilterUtils.areProperSubset(pf1, pf1));
      assertTrue("are subset: " + pf2 + ", " + pf2, FilterUtils.areProperSubset(pf2, pf2));
      assertTrue("are subset: " + pf3 + ", " + pf3, FilterUtils.areProperSubset(pf3, pf3));

      assertTrue("not disjoint: " + pf1 + ", " + pf1, !FilterUtils.areDisjoint(pf1, pf1));
      assertTrue("not disjoint: " + pf2 + ", " + pf2, !FilterUtils.areDisjoint(pf2, pf2));
      assertTrue("not disjoint: " + pf3 + ", " + pf3, !FilterUtils.areDisjoint(pf3, pf3));
    }
    
    public void testByProperty() throws Exception {
      assertTrue("are disjoint: " + pf4 + ", " + pf5, FilterUtils.areDisjoint(pf4, pf5));
      assertTrue("not disjoint: " + pf4 + ", " + pf6, !FilterUtils.areDisjoint(pf4, pf6));
      assertTrue("not disjoint: " + pf5 + ", " + pf6, !FilterUtils.areDisjoint(pf5, pf6));
      
      assertTrue("not subset: " + pf4 + ", " + pf5, !FilterUtils.areProperSubset(pf4, pf5));
      assertTrue("not subset: " + pf4 + ", " + pf6, !FilterUtils.areProperSubset(pf4, pf6));
      assertTrue("not subset: " + pf5 + ", " + pf6, !FilterUtils.areProperSubset(pf5, pf6));

      assertTrue("are subset: " + pf4 + ", " + pf4, FilterUtils.areProperSubset(pf4, pf4));
    }
    
    public void testHasByProperty() throws Exception {
      assertTrue("not disjoint: " + pf1 + ", " + pf4, !FilterUtils.areDisjoint(pf1, pf4));
      assertTrue("not disjoint: " + pf2 + ", " + pf4, !FilterUtils.areDisjoint(pf2, pf4));
      assertTrue("not disjoint: " + pf3 + ", " + pf4, !FilterUtils.areDisjoint(pf3, pf4));
      
      assertTrue("not subset: " + pf1 + ", " + pf4, !FilterUtils.areProperSubset(pf1, pf4));
      assertTrue("are subset: " + pf4 + ", " + pf1, FilterUtils.areProperSubset(pf4, pf1));
      assertTrue("are subset: " + pf5 + ", " + pf1, FilterUtils.areProperSubset(pf5, pf1));
      assertTrue("not subset: " + pf6 + ", " + pf1, !FilterUtils.areProperSubset(pf6, pf1));

      assertTrue("not subset: " + pf2 + ", " + pf4, !FilterUtils.areProperSubset(pf2, pf4));
      assertTrue("not subset: " + pf4 + ", " + pf2, !FilterUtils.areProperSubset(pf4, pf2));
      assertTrue("not subset: " + pf5 + ", " + pf2, !FilterUtils.areProperSubset(pf5, pf2));
      assertTrue("are subset: " + pf6 + ", " + pf2, FilterUtils.areProperSubset(pf6, pf2));
    }
    
    public void testAndOptimizeAllNone() {
      optimizeExact(all_and_all, FilterUtils.all());
      optimizeExact(FilterUtils.and(all_and_all, all_and_all), FilterUtils.all());
      
      optimizeExact(all_and_none, FilterUtils.none());
      optimizeExact(FilterUtils.and(all_and_none, FilterUtils.all()), FilterUtils.none());
      
      optimizeExact(none_and_none, FilterUtils.none());
    }
    
    public void testOrOptimizeAllNone() {
      optimizeExact(all_or_all, FilterUtils.all());
      optimizeExact(FilterUtils.or(all_or_all, all_or_all), FilterUtils.all());
      
      optimizeExact(all_or_none, FilterUtils.all());
      optimizeExact(FilterUtils.or(all_or_none, FilterUtils.none()), FilterUtils.all());
      
      optimizeExact(none_or_none, FilterUtils.none());
    }
    
    public void testAndOrAllNone() {
      optimizeExact(FilterUtils.or(all_and_all, FilterUtils.all()), FilterUtils.all());
      optimizeExact(FilterUtils.or(all_and_none, FilterUtils.all()), FilterUtils.all());

      optimizeExact(FilterUtils.or(all_and_all, FilterUtils.none()), FilterUtils.all());
      optimizeExact(FilterUtils.or(all_and_none, FilterUtils.none()), FilterUtils.none());

      optimizeExact(FilterUtils.and(all_or_all, FilterUtils.all()), FilterUtils.all());
    }
    
    public void testUseCases() {
      // schemas for our database features
      FeatureFilter transcript = FilterUtils.byType("transcript");
      FeatureFilter exon = FilterUtils.byType("exon");
      FeatureFilter repeat = FilterUtils.byType("repeat");
      
      AnnotationType.Impl tsType = new AnnotationType.Impl();
      tsType.setConstraints("transcript.id", new PropertyConstraint.ByClass(String.class), CardinalityConstraint.ONE);
      tsType.setDefaultConstraints(PropertyConstraint.NONE, CardinalityConstraint.NONE);
      FeatureFilter tsID = FilterUtils.byAnnotationType(tsType);
      
      AnnotationType.Impl exType = new AnnotationType.Impl();
      exType.setConstraints("id", new PropertyConstraint.ByClass(String.class), CardinalityConstraint.ONE);
      exType.setDefaultConstraints(PropertyConstraint.NONE, CardinalityConstraint.NONE);
      FeatureFilter exID = FilterUtils.byAnnotationType(exType);
      
      AnnotationType.Impl rpType = new AnnotationType.Impl();
      rpType.setConstraints("id", PropertyConstraint.ANY, CardinalityConstraint.ONE);
      rpType.setDefaultConstraints(PropertyConstraint.NONE, CardinalityConstraint.NONE);
      FeatureFilter repeatID = FilterUtils.byAnnotationType(rpType);
      
      FeatureFilter tsSchema = FilterUtils.and(transcript, tsID);
      FeatureFilter exSchema = FilterUtils.and(exon, exID);
      FeatureFilter reSchema = FilterUtils.and(repeat, repeatID);
      
      FeatureFilter dbFilter = FilterUtils.or(
        new FeatureFilter[] { tsSchema, exSchema, reSchema }
      );
      
      // pull out a feature by transcript.id
      FeatureFilter aTranscript = FilterUtils.byAnnotation("transcript.id", "ts:42");
      
      // let the fun commence
      optimizeEquals(FilterUtils.and(tsSchema, aTranscript), FilterUtils.and(transcript, aTranscript));
      optimizeExact(FilterUtils.and(exSchema, aTranscript), FilterUtils.none());
      optimizeExact(FilterUtils.and(reSchema, aTranscript), FilterUtils.none());
      optimizeEquals(dbFilter, dbFilter);
      optimizeEquals(FilterUtils.and(dbFilter, aTranscript), FilterUtils.and(transcript, aTranscript));
    }
    
    public void testAndOrMethods() {
      checkEquals(
        FilterUtils.and(FilterUtils.and(tf1, tf2), tf3),
        FilterUtils.and(new FeatureFilter[] { tf1, tf2, tf3 })
      );

      checkEquals(
        FilterUtils.and(tf1, FilterUtils.and(tf2, tf3)),
        FilterUtils.and(new FeatureFilter[] { tf1, tf2, tf3 })
      );
      
      checkEquals(
        FilterUtils.or(FilterUtils.or(tf1, tf2), tf3),
        FilterUtils.or(new FeatureFilter[] { tf1, tf2, tf3 })
      );

      checkEquals(
        FilterUtils.or(tf1, FilterUtils.or(tf2, tf3)),
        FilterUtils.or(new FeatureFilter[] { tf1, tf2, tf3 })
      );
    }
    
    public void testAncestorDescendants() {
        FeatureFilter type_foo = new FeatureFilter.ByType("foo");
        FeatureFilter type_bar = new FeatureFilter.ByType("bar");
        FeatureFilter test = new FeatureFilter.ByAncestor(new FeatureFilter.OnlyDescendants(type_foo));
        
        assertTrue(
            FilterUtils.areDisjoint(
                type_bar,
                test
            )
        );
        assertTrue(
            !FilterUtils.areDisjoint(
                type_foo,
                test
            )
        );
    }
    
    public void testOrAncestorDescendants() {
        FeatureFilter type_foo = new FeatureFilter.ByType("foo");
        FeatureFilter type_bar = new FeatureFilter.ByType("bar");
        FeatureFilter test = new FeatureFilter.ByAncestor(new FeatureFilter.OnlyDescendants(type_foo));
        
        assertTrue(
            FilterUtils.areDisjoint(
                type_bar,
                test
            )
        );
        assertTrue(
            !FilterUtils.areDisjoint(
                type_foo,
                test
            )
        );
    }
    
    public void testAncestorLeafChildren() {
        FeatureFilter type_foo = new FeatureFilter.ByType("foo");
        FeatureFilter type_bar = new FeatureFilter.ByType("bar");
        FeatureFilter test = new FeatureFilter.ByAncestor(
                new FeatureFilter.OnlyChildren(
                        new FeatureFilter.And(
                                type_foo,
                                FeatureFilter.leaf
                        )
                )
        );
        FeatureFilter test_nonleaf = new FeatureFilter.ByAncestor(
                new FeatureFilter.OnlyChildren(
                                type_foo
                )
        );
        
        assertTrue(
            FilterUtils.areDisjoint(
                type_bar,
                test
            )
        );
        assertTrue(
            !FilterUtils.areDisjoint(
                type_bar,
                test_nonleaf
            )
        );
        assertTrue(
            !FilterUtils.areDisjoint(
                type_foo,
                test
            )
        );
    }
    
    public void testAncestorComplexChildren() {
        FeatureFilter type_foo = new FeatureFilter.ByType("foo");
        FeatureFilter type_bar = new FeatureFilter.ByType("bar");
        FeatureFilter type_baz = new FeatureFilter.ByType("baz");
        FeatureFilter test = new FeatureFilter.ByAncestor(
                new FeatureFilter.OnlyChildren(
                        new FeatureFilter.And(
                                type_foo,
                                new FeatureFilter.OnlyChildren(
                                        new FeatureFilter.And(
                                                type_baz,
                                                FeatureFilter.leaf
                                        )
                                )
                        )
                )
        );
        
        assertTrue(
            FilterUtils.areDisjoint(
                type_bar,
                test
            )
        );
        assertTrue(
            !FilterUtils.areDisjoint(
                type_foo,
                test
            )
        );
        assertTrue(
            !FilterUtils.areDisjoint(
                type_baz,
                test
            )
        );
    }
    
    public void testParentChildren() {
        FeatureFilter type_foo = new FeatureFilter.ByType("foo");
        FeatureFilter type_bar = new FeatureFilter.ByType("bar");
        FeatureFilter type_baz = new FeatureFilter.ByType("baz");
        FeatureFilter test = new FeatureFilter.ByParent(
                new FeatureFilter.OnlyChildren(
                        new FeatureFilter.And(
                                type_foo,
                                new FeatureFilter.OnlyChildren(
                                        new FeatureFilter.And(
                                                type_baz,
                                                FeatureFilter.leaf
                                        )
                                )
                        )
                )
        );
        
        assertTrue(
            FilterUtils.areDisjoint(
                type_bar,
                test
            )
        );
        assertTrue(
            !FilterUtils.areDisjoint(
                type_foo,
                test
            )
        );
        assertTrue(
            FilterUtils.areDisjoint(
                type_baz,
                test
            )
        );
    }
    
    public void testDisjointAnnotationTypes() {
        AnnotationType atype = new AnnotationType.Impl();
        atype.setDefaultConstraints(PropertyConstraint.NONE, CardinalityConstraint.ZERO);
        atype.setConstraints(
            "foo",
            PropertyConstraint.ANY,
            CardinalityConstraint.ONE
        );
        assertTrue(FilterUtils.areDisjoint(
            new FeatureFilter.ByAnnotationType(atype),
            new FeatureFilter.ByAnnotation("bar", "some_value")
        ));
        assertTrue(!FilterUtils.areDisjoint(
            new FeatureFilter.ByAnnotationType(atype),
            new FeatureFilter.ByAnnotation("foo", "some_value")
        ));
    }
    
    public void testDisjointAnnotationTypesContains() {
        AnnotationType atype = new AnnotationType.Impl();
        atype.setDefaultConstraints(PropertyConstraint.NONE, CardinalityConstraint.ZERO);
        atype.setConstraints(
            "foo",
            PropertyConstraint.ANY,
            CardinalityConstraint.ONE
        );
        assertTrue(FilterUtils.areDisjoint(
            new FeatureFilter.ByAnnotationType(atype),
            new FeatureFilter.AnnotationContains("bar", "some_value")
        ));
        assertTrue(!FilterUtils.areDisjoint(
            new FeatureFilter.ByAnnotationType(atype),
            new FeatureFilter.AnnotationContains("foo", "some_value")
        ));
    }
    
    private void optimizeExact(FeatureFilter raw, FeatureFilter target) {
      FeatureFilter result = FilterUtils.optimize(raw);
      assertTrue("optimize: " + raw + " should be " + target + " but is " + result, result == target);
    }
    
    private void optimizeEquals(FeatureFilter raw, FeatureFilter target) {
      FeatureFilter result = FilterUtils.optimize(raw);
      assertTrue("optimize: " + raw + " should be " + target + " but is " + result, result.equals(target));
    }
    
    private void checkEquals(FeatureFilter a, FeatureFilter b) {
      assertTrue("equal: " + a + ", " + b, a.equals(b));
    }
}
