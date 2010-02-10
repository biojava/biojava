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

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * StAX handler for XFF strandedFeature type.  By default, XFFFeatureSetHandler
 * uses this handler for all feature elements which have the <code>strand</code>
 * attribute.
 *
 * <p>
 * Like the basic <code>FeatureHandler</code>, this class can be subclassed
 * to give handlers for more specialized feature types.
 * </p>
 *
 * @author Thomas Down
 * @since 1.2
 */

public class StrandedFeatureHandler extends FeatureHandler {
    boolean inFeature = false;

    public static final XFFPartHandlerFactory STRANDEDFEATURE_HANDLER_FACTORY = new XFFPartHandlerFactory() {
	    public StAXContentHandler getPartHandler(XFFFeatureSetHandler xffenv) {
		return new StrandedFeatureHandler(xffenv);
	    }
	} ;

    public StrandedFeatureHandler(XFFFeatureSetHandler xffenv) {
	super(xffenv);
    }

    protected Feature.Template createFeatureTemplate() {
	return new StrandedFeature.Template();
    }

    protected StrandedFeature.Template getStrandedFeatureTemplate() {
	return (StrandedFeature.Template) getFeatureTemplate();
    }

    public void startElement(String nsURI,
			     String localName,
			     String qName,
			     Attributes attrs,
			     DelegationManager dm)
	 throws SAXException
    {
	if (!inFeature) {
	    // This is the root startElement.  Check the strand= attribute.

	    String strands = attrs.getValue("strand");
	    StrandedFeature.Template ft = getStrandedFeatureTemplate();
	    if (strands != null) {
		if (strands.equals("+")) {
		    ft.strand = StrandedFeature.POSITIVE;
		} else if (strands.equals("-")) {
		    ft.strand = StrandedFeature.NEGATIVE;
		} else {
		    ft.strand = StrandedFeature.UNKNOWN;
		}
	    } else {
		ft.strand = StrandedFeature.UNKNOWN;
	    }

	    inFeature = true;
	}

	// Pass everything on to the basic feature parser.

	super.startElement(nsURI, localName, qName, attrs, dm);
    }
}
