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

import org.biojava.bio.seq.io.ParseException;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * StAX handler for xff:prop detail elements.  prop details are simple
 * tag-value text, and they are added directly as properties of the feature.
 *
 * @author Thomas Down
 * @since 1.2
 */

public class PropDetailHandler extends StringElementHandlerBase {
    public static final XFFPartHandlerFactory PROPDETAIL_HANDLER_FACTORY = new XFFPartHandlerFactory() {
	    public StAXContentHandler getPartHandler(XFFFeatureSetHandler xffenv) {
		return new PropDetailHandler(xffenv);
	    }
	} ;

    private String key;
    private XFFFeatureSetHandler xffenv;

    public PropDetailHandler(XFFFeatureSetHandler xffenv) {
	super();
	this.xffenv = xffenv;
    }

    public void startElement(String nsURI,
			     String localName,
			     String qName,
			     Attributes attrs,
			     DelegationManager dm)
	 throws SAXException
    {
	key = attrs.getValue("key");
	if (key == null) {
	    throw new SAXException("No key attribute for xff:prop");
	}
	super.startElement(nsURI, localName, qName, attrs, dm);
    }

    protected void setStringValue(String s)
        throws SAXException
    {
	try {
	    xffenv.getFeatureListener().addFeatureProperty(key, s);
	} catch (ParseException pex) {
	    throw new SAXException(pex);
	}
    }
}
