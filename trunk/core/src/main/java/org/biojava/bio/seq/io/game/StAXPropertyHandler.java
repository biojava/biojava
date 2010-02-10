/**
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

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * StAX handler shamelessly ripped off from Thomas Down's
 * XFFFeatureSetHandler.
 *
 * <strong>NOTE</strong> This class is not thread-safe -- it
 * must only be used for one parse at any time.
 *
 * This class is the basis for classes that do not create
 * a new feature but modify an existing feature.
 *
 * It is not compulsory for property handlers to subclass
 * this class but those that don't but wish to use the
 * handler stack facility need to use the StaxFeatureHandler's
 * push and pop methods.
 *
 * @author Thomas Down
 * @author David Huen
 *
 * @since 1.8
 */
public class StAXPropertyHandler extends StAXContentHandlerBase
{
  private String myLocalName;
  private boolean hasCallback = false;
  private boolean inElement = false;
  private boolean setOnceFired = false;

  // class variables
  protected SeqIOListener featureListener;
  private List handlers;
  protected StAXFeatureHandler staxenv;
  private int baseLevel = 0;

  {
    handlers = new ArrayList();
  }

  // there should be a factory method here to make this class

  // constructor
  // because every StAXPropertyHandler was ultimately
  // invoked from a StAXFeatureHandler via delegation
  // the staxenv will point at the StAXFeatureHandler
  // at the base of the current chain of StAXPropertyHandler,
  // which is not necessarily the first root element
  // StAXFeatureHandler.
  StAXPropertyHandler(StAXFeatureHandler staxenv) {
    // cache environmnet
    this.staxenv = staxenv;
  }

/**
 * Sets the element name that the class responds to.
 */
  public void setHandlerCharacteristics(String localName, boolean hasCallback) 
  {    
  if (!setOnceFired) {
      myLocalName = localName;
      this.hasCallback = hasCallback;
      setOnceFired = true;
    }
    else
      System.err.println("setHandlerChracteristics called twice on same handler");  
  }


/**
 * get iterator for current stack starting at the position
 * below mine.
 */
  protected ListIterator getHandlerStackIterator() 
      throws ParseException {
    if (baseLevel >= 1)
      return staxenv.getHandlerStackIterator(baseLevel-1);
    else
      throw new ParseException("getHandlerStackIterator while at bottom of stack.");
  }

  // Class to implement bindings
  class Binding {
    final ElementRecognizer recognizer;
    final StAXHandlerFactory handlerFactory;
    Binding(ElementRecognizer er, StAXHandlerFactory hf)
    {
      recognizer = er;
      handlerFactory = hf;
    }
  }

  // method to add a handler
  // we do not distinguish whither it is a feature or property
  // handler.  The factory method creates the right type subclassed
  // from the correct type of handler
  protected void addHandler(
                   ElementRecognizer rec, 
                   StAXHandlerFactory handler)
  {
//    System.out.println("StAXPropertyHandler.addHandler called.");
    handlers.add(new Binding(rec, handler));
//    System.out.println("StAXPropertyHandler.addHandler left.");
//    System.out.println(" ");
  }

/**
 * Element-specific handler.
 * Subclass this to do something useful!
 */
  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
         throws SAXException
  {
  }

/**
 * Override this to do any processing required but call this
 * prior to returning.  Delegation occurs here!
 *
 */
  public void startElement(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs,
                DelegationManager dm)
	 throws SAXException
  {
//    System.out.println("StAXPropertyHandler.startElement localName: " + localName);
    // perform delegation
    for (int i = handlers.size() - 1; i >= 0; --i) {
      Binding b = (Binding) handlers.get(i);
      if (b.recognizer.filterStartElement(nsURI, localName, qName, attrs)) {
      dm.delegate(b.handlerFactory.getHandler(staxenv));
      return;
      }
    }

    // is this for me?
    if (!(myLocalName.equals(localName)) ) return;

    if (!inElement) {
      // save current stack position just in case I want to search downwards.
      baseLevel = staxenv.getLevel();

      if (hasCallback) staxenv.push(this);
      
      inElement = true;
    }

    if (inElement) startElementHandler(nsURI, localName, qName, attrs);
  }

/**
 * Element specific exit handler
 * Subclass to do anything useful.
 */
  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
              throws SAXException
  {
  }

  public void endElement(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
              throws SAXException
  {
    // is this mine?
    if (!(myLocalName.equals(localName)) ) return;

//    System.out.println("StAXPropertyHandler.endElement localName: " + localName);
    // do the necessary before exit
    if (inElement) {
      // element specific handling
      endElementHandler(nsURI, localName, qName, handler);

      if (hasCallback)
        if (setOnceFired) {
          staxenv.pop();
          setOnceFired = false;
        }
   
      inElement = false;
    }
  }
}
