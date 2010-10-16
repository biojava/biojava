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

/**
 * Tests for SimpleAssembly.  By dependancy, this also
 * tests ProjectedFeatureHolder and SimpleAssembly.
 *
 * @author Thomas Down
 * @since 1.3
 */

public class MergeFeatureHolderTest extends TestCase
{
    public MergeFeatureHolderTest(String name) {
        super(name);
    }
    
    public void testTopLevelAccess() 
        throws Exception
    {
        MergeFeatureHolder mfh = new MergeFeatureHolder();
        mfh.addFeatureHolder(new SimpleFeatureHolder(
                new FeatureFilter.And(
                        FeatureFilter.top_level,
                        new FeatureFilter.And(
                                new FeatureFilter.ByType("foo"),
                                FeatureFilter.leaf
                        )
               )
        ));
        mfh.addFeatureHolder(new FailOnAccessFeatureHolder(
                new FeatureFilter.And(
                        FeatureFilter.top_level,
                        new FeatureFilter.And(
                                new FeatureFilter.ByType("bar"),
                                new FeatureFilter.OnlyChildren(
                                        new FeatureFilter.And(
                                                new FeatureFilter.ByType("foo"),
                                                FeatureFilter.leaf
                                        )
                                )
                        )
               )
        ));
        
        mfh.filter(new FeatureFilter.ByType("foo"), false).countFeatures();
    }
    
    public void testRecursiveAccess() 
        throws Exception
    {
        MergeFeatureHolder mfh = new MergeFeatureHolder();
        mfh.addFeatureHolder(new SimpleFeatureHolder(
                new FeatureFilter.And(
                        FeatureFilter.top_level,
                        new FeatureFilter.And(
                                new FeatureFilter.ByType("foo"),
                                FeatureFilter.leaf
                        )
               )
        ));
        mfh.addFeatureHolder(new FailOnAccessFeatureHolder(
                new FeatureFilter.And(
                        FeatureFilter.top_level,
                        new FeatureFilter.And(
                                new FeatureFilter.ByType("bar"),
                                new FeatureFilter.OnlyChildren(
                                        new FeatureFilter.And(
                                                new FeatureFilter.ByType("baz"),
                                                FeatureFilter.leaf
                                        )
                                )
                        )
               )
        ));
        
        mfh.filter(new FeatureFilter.ByType("quux")).countFeatures();
    }
    
    private class FailOnAccessFeatureHolder extends LazyFeatureHolder {
        public FailOnAccessFeatureHolder(FeatureFilter schema) {
            super(schema);
        }
        
        protected FeatureHolder createFeatureHolder() {
            throw new RuntimeException("This FeatureHolder should not have been accessed");
        }
    }
}
