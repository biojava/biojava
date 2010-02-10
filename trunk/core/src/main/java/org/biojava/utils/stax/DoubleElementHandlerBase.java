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
 * Copyright for this code is held jodoublely by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package org.biojava.utils.stax;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * <p>
 * StAX handler for any element which just contains a string representation of
 * a double.
 * </p>
 *
 * <p>
 * This calss collects the string data, and when it is complete, passes it to
 * the (abstract) setDoubleValue method.  Typical use of this class is as
 * a base for a small (often anonymous) class which takes the double value
 * and stores it in some variable.
 * </p>
 *
 * @author Matthew Pocock
 * @author Greg Cox
 * @since 1.2
 */

public abstract class DoubleElementHandlerBase extends StAXContentHandlerBase {
  private int level = 0;
  private StringBuffer data = new StringBuffer();

  public void startElement(
    String nsURI,
    String localName,
    String qName,
    Attributes attrs,
    DelegationManager dm
  ) throws SAXException {
    level++;
    if (level > 1) {
      throw new SAXException("Found child element when expecting character data");
    }
  }

  public void endElement(
    String nsURI,
    String localName,
    String qName,
    StAXContentHandler handler
  ) throws SAXException {
    level--;
    if (level == 0) {
      try {
        setDoubleValue(Double.parseDouble((data.substring(0)).trim()));
      } catch (NumberFormatException nfe) {
        throw new SAXException(nfe);
      }
    }
  }

  public void characters(char[] ch, int start, int end) throws SAXException {
    data.append(ch, start, end);
  }

  /**
   * <p>
   * Override this method to do something useful with the
   * double we collect.
   * </p>
   *
   * <p>
   * This method will be invoked by endElement with the fully parsed double.
   * </p>
   *
   * @param val the fully parsed double
   * @throws SAXException if for any reason the double is not palatable
   */
  protected abstract void setDoubleValue(double val) throws SAXException;
}
