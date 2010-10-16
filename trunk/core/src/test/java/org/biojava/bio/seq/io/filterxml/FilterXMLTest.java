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

package org.biojava.bio.seq.io.filterxml;

import java.io.CharArrayReader;
import java.io.CharArrayWriter;
import java.io.PrintWriter;

import javax.xml.parsers.SAXParserFactory;

import junit.framework.TestCase;

import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.FramedFeature;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.stax.SAX2StAXAdaptor;
import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;
import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;

/**
 * Test round-tripping of FeatureFilters via the FilterXML language 
 *
 * @author Thomas Down
 */
public class FilterXMLTest extends TestCase
{
    XMLFilterWriter filterWriter;
    XMLFilterHandler filterHandler;
    
    public FilterXMLTest(String theString)
    {
        super(theString);
        this.filterWriter = new XMLFilterWriter();
        this.filterHandler = new XMLFilterHandler();
    }

    private FeatureFilter roundTripFilter(FeatureFilter ff) 
        throws Exception
    {
        CharArrayWriter caw = new CharArrayWriter();
        PrintWriter pw = new PrintWriter(caw);
        XMLWriter xw = new PrettyXMLWriter(pw);
        filterWriter.writeFilter(ff, xw);
        pw.flush();
        
        CharArrayReader car = new CharArrayReader(caw.toCharArray());
        SAXParserFactory spf = SAXParserFactory.newInstance();
	    spf.setNamespaceAware(true);
	    XMLReader parser = spf.newSAXParser().getXMLReader();
        XMLFilterHandler.FilterHandler handler = filterHandler.getStAXContentHandler();
        parser.setContentHandler(new SAX2StAXAdaptor(handler));
        parser.parse(new InputSource(car));
        
        FeatureFilter rtff = handler.getFeatureFilter();
        assertTrue("Testing that " + ff.toString() + " contains " + rtff.toString(), FilterUtils.areProperSubset(rtff, ff));
        assertTrue("Testing that " + rtff.toString() + " contains " + ff.toString(), FilterUtils.areProperSubset(ff, rtff));
        return rtff;
    }
    
    public void testByType()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.ByType("foo"));
    }
    
    public void testNot()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.Not(new FeatureFilter.ByType("bar")));
    }
    
    public void testAnd()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.And(new FeatureFilter.ByType("foo"), new FeatureFilter.BySource("bar")));
    }
    
    public void testMultiAnd()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.And(new FeatureFilter.And(new FeatureFilter.ByType("foo"), new FeatureFilter.BySource("bar")), new FeatureFilter.ByComponentName("AL121903")));
        roundTripFilter(new FeatureFilter.And(new FeatureFilter.ByType("foo"), new FeatureFilter.And(new FeatureFilter.BySource("bar"), new FeatureFilter.ByComponentName("AL121903"))));
    }
    
    public void testLocation()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.OverlapsLocation(new RangeLocation(500, 2000)));
        Location complex = LocationTools.union(
            new RangeLocation(100, 200),
            new RangeLocation(500, 600)
        );
        roundTripFilter(new FeatureFilter.ContainedByLocation(complex));
    }
    
    public void testHasAnnotation()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.HasAnnotation("foo"));
    }
    
    public void testByAnnotation()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.ByAnnotation("foo", "bar"));
    }
    
    public void testStrandFilter()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.StrandFilter(StrandedFeature.POSITIVE));
        roundTripFilter(new FeatureFilter.StrandFilter(StrandedFeature.NEGATIVE));
        roundTripFilter(new FeatureFilter.StrandFilter(StrandedFeature.UNKNOWN));
    }
    
    public void testFrameFilter()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.FrameFilter(FramedFeature.FRAME_1));
    }
    
    public void testByClass()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.ByClass(ComponentFeature.class));
    }
    
    public void testByParent()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.ByParent(new FeatureFilter.ByType("foo")));
    }
    
    public void testAnnotationContains()
        throws Exception
    {
        roundTripFilter(new FeatureFilter.AnnotationContains("foo.bar", "baz"));
    }
}
