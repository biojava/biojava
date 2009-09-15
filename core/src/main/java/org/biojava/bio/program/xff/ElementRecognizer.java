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

import org.xml.sax.Attributes;

/**
 * Simple interface for filtering SAX/StAX startElement events.
 *
 * <p>
 * A number of standard implementations are provided for your convenience. To
 * implement your own filters, just implement the filterStartElement method.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.2
 */

public interface ElementRecognizer {
  /**
   * Recognize an element based upon the start element parameters.
   *
   * @param nsURI     the uri of the element to filter
   * @param localName the local name of the element to filter
   * @param qName     the qName of the element to filter
   * @param attrs     the attributes associated with the element to filter
   * @return          true if this element is accepted, false otherwise
   */
  public boolean filterStartElement(String nsURI, String localName, String qName, Attributes attrs);

  public static final ElementRecognizer ALL = new AllElementRecognizer();

  public static class AllElementRecognizer implements ElementRecognizer {
    public boolean filterStartElement(String nsURI,
                                      String localName,
                                      String qName,
                                      Attributes attrs)
    {
      return true;
    }

    public String toString()
    {
      return "ALL";
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

    public String toString()
    {
      return "HasAttribute[name=" + attributeName + "]";
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

    public String toString()
    {
      return "ByNSName[uri=" + nsURI + " name=" + localName + "]";
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

    public String toString()
    {
      return "ByLocalName[name=" + localName + "]";
    }
  }
}
