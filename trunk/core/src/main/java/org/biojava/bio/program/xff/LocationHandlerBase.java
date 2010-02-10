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

import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * Handler to the XFF location type.  To use this, write a simple
 * subclass which takes the location and stores it somewhere sensible.
 *
 * @author Thomas Down
 * @since 1.2
 */

public abstract class LocationHandlerBase extends StAXContentHandlerBase {
    private int level = 0;
    private Location location = Location.empty;

    public void startElement(String nsURI,
			     String localName,
			     String qName,
			     Attributes attrs,
			     DelegationManager dm)
	 throws SAXException
    {
	level++;
	if (level == 0) {
	    // Top-level.  Ignore
	} else {
	    if (localName.equals("span")) {
		try {
		    String starts = attrs.getValue("start");
		    String stops = attrs.getValue("stop");
		    
		    if (starts == null || stops == null) {
			throw new SAXException("Missing start/stop attribute in location span");
		    }
		    
		    Location span = new RangeLocation(Integer.parseInt(starts),
						      Integer.parseInt(stops));
		    location = location.union(span);
		} catch (NumberFormatException ex) {
		    throw new SAXException("Error parsing location", ex);
		}
	    }
	}
    }

    public void endElement(String nsURI,
			   String localName,
			   String qName,
			   StAXContentHandler handler)
	throws SAXException
    {
	level--;
	if (level == 0) {
	    setLocationValue(location);
	}
    }

    /**
     * Override this method to do something useful with the
     * location we collect.  Maybe we should do this by delegation
     * rather than extension.
     */

    protected abstract void setLocationValue(Location l) throws SAXException;
}
