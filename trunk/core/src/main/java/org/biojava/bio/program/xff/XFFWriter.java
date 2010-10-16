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

package org.biojava.bio.program.xff;

import java.io.IOException;
import java.util.Iterator;

import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.xml.XMLWriter;

/**
 * Simple XFF writing code, ripped off from Dazzle 0.08.  It probably
 * needs a re-write to allow custom types to be used properly, but
 * it has enough hooks for now...
 *
 * @author Thomas Down
 */

public class XFFWriter {
    private XFFHelper helper;

    public XFFWriter() {
	helper = new BasicXFFHelper();
    }

    public XFFWriter(XFFHelper h) {
	helper = h;
    }

    private void emitFeature_xff(Feature f,
				 XMLWriter xw)
	throws IOException
    {
	if (f instanceof ComponentFeature) {
	    ComponentFeature cf = (ComponentFeature) f;

	    String strand = (cf.getStrand() == StrandedFeature.POSITIVE ? "+" : "-");
	    xw.openTag("componentFeature");
	    xw.attribute("strand", strand);
	    String id = helper.getFeatureID(cf);
	    if (id != null) {
		xw.attribute("id", id);
	    }
	    emitFeature_xff_standardBits(f, xw);
	    
	    xw.openTag("componentID");
	    xw.print(cf.getComponentSequence().getName());
	    xw.closeTag("componentID");
	    
	    xw.openTag("componentLocation");
	    emitLocation_xff(cf.getComponentLocation(), xw);
	    xw.closeTag("componentLocation");

	    emitDetails(xw, cf);

	    xw.closeTag("componentFeature");
	} else if (f instanceof StrandedFeature && ((StrandedFeature) f).getStrand() != StrandedFeature.UNKNOWN) {
	    StrandedFeature sf = (StrandedFeature) f;

	    String strand = (sf.getStrand() == StrandedFeature.POSITIVE ? "+" : "-");
	    xw.openTag("strandedFeature");
	    xw.attribute("strand", strand);
	    String id = helper.getFeatureID(sf);
	    if (id != null) {
		xw.attribute("id", id);
	    }
	    emitFeature_xff_standardBits(f, xw);
	    if (f.countFeatures() > 0)
		emitFeature_xff_featureSet(f, xw, false);

	    emitDetails(xw, f);

	    xw.closeTag("strandedFeature");
	} else {
	    xw.openTag("feature");
	    String id = helper.getFeatureID(f);
	    if (id != null) {
		xw.attribute("id", id);
	    }
	    emitFeature_xff_standardBits(f, xw);
	    if (f.countFeatures() > 0)
		emitFeature_xff_featureSet(f, xw, false);
	    
	    emitDetails(xw, f);

	    xw.closeTag("feature");
	}
    }

    private void emitDetails(XMLWriter xw,
			     Feature f)
	throws IOException
    {
	xw.openTag("details");
	helper.writeDetails(xw, f);
	xw.closeTag("details");
    }

    private void emitFeature_xff_standardBits(Feature f,
					      XMLWriter xw)
	throws IOException
    {
	xw.openTag("type");
	xw.print(f.getType());
	xw.closeTag("type");

	xw.openTag("source");
	xw.print(f.getSource());
	xw.closeTag("source");

	xw.openTag("location");
	emitLocation_xff(f.getLocation(), xw);
	xw.closeTag("location");
    }

    private void emitLocation_xff(Location l,
				  XMLWriter xw)
	throws IOException
    {
	for (Iterator i = l.blockIterator(); i.hasNext(); ) {
	    Location block = (Location) i.next();
	    xw.openTag("span");
	    xw.attribute("start", "" + block.getMin());
	    xw.attribute("stop", "" + block.getMax());
	    xw.closeTag("span");
	}
    }

    
    private void emitFeature_xff_featureSet(FeatureHolder fh,
					    XMLWriter xw,
					    boolean explicitNamespace)
	throws IOException
    {
	xw.openTag("featureSet");
	if (explicitNamespace) {
		xw.declareNamespace(XFFTools.XFF_NS, "");
		xw.declareNamespace(XFFTools.XFF_BIOJAVA_NS, "biojava");
	}
	
	for (Iterator fi = fh.features(); fi.hasNext(); ) {
	    Feature f = (Feature) fi.next();
	    emitFeature_xff(f, xw);
	}
	xw.closeTag("featureSet");
    }

    public void writeFeatureSet(FeatureHolder fh,
				XMLWriter xw)
	throws IOException
    {
	emitFeature_xff_featureSet(fh, xw, true);
    }
}
