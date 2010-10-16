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

package org.biojava.bio.seq.io.game;

import org.xml.sax.Attributes;

/**
 * Simple interface for filtering SAX/StAX startElement events
 *
 * @author Thomas Down
 * @since 1.2
 */

// this class really should be refactored elsewhere.

public interface ElementRecognizer {
    public static final ElementRecognizer ALL = new AllElementRecognizer();

    public static class AllElementRecognizer implements ElementRecognizer {
	public boolean filterStartElement(String nsURI,
					      String localName,
					      String qName,
					      Attributes attrs)
	{
	    return true;
	}
    }

    /**
     * Filter elements on the existence of a specified attribute.
     */

    public static class HasAttribute implements ElementRecognizer {
	private String attributeName;

	public HasAttribute(String name) {
	    attributeName = name;
	}

	public boolean filterStartElement(String nsURI,
					  String localName,
					  String qName,
					  Attributes attrs)
	{
	    return (attrs.getValue(attributeName) != null);
	}
    }

    /**
     * Filter elements by name and namespace.
     */

    public static class ByNSName implements ElementRecognizer {
	private String nsURI;
	private String localName;

	public ByNSName(String nsURI, String localName) {
	    this.nsURI = nsURI;
	    this.localName = localName;
	}

	public boolean filterStartElement(String nsURI,
					  String localName,
					  String qName,
					  Attributes attrs)
	{
	    return (localName.equals(this.localName) && nsURI.equals(this.nsURI));
	}
    }

    /**
     * Filter elements by local name (not recommended).
     */

    public static class ByLocalName implements ElementRecognizer {
	private String localName;

	public ByLocalName(String localName) {
	    this.localName = localName;
	}

	public boolean filterStartElement(String nsURI,
					  String localName,
					  String qName,
					  Attributes attrs)
	{
	    return (localName.equals(this.localName));
	}
    }

    public boolean filterStartElement(String nsURI, String localName, String qName, Attributes attrs);
}
