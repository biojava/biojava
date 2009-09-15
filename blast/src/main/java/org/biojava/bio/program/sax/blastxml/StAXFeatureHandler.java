/*
 *  BioJava development code This code may be freely distributed and modified
 *  under the terms of the GNU Lesser General Public Licence. This should be
 *  distributed with the code. If you do not have a copy, see:
 *  http://www.gnu.org/copyleft/lesser.html Copyright for this code is held
 *  jointly by the individual authors. These should be listed in
 *
 *@author    doc comments. For more information on the BioJava project and its
 *      aims, or to join the biojava-l mailing list, visit the home page at:
 *      http://www.biojava.org/
 */

package org.biojava.bio.program.sax.blastxml;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;
import org.xml.sax.SAXException;

/**
 *  StAX handler shamelessly ripped off from Thomas Down's XFFFeatureSetHandler.
 *  It was modified for greater generality. <strong>NOTE</strong> This class is
 *  not thread-safe -- it must only be used for one parse at any time.
 *
 *@author     Thomas Down
 *@author     David Huen
 *@created    19 January 2002
 *@since      1.8
 */

class StAXFeatureHandler extends StAXContentHandlerBase {
    // class variables
    /**
     * declare namespace
     */
    static final String biojavaUri = "http://www.biojava.org";
//    static final String biojavaUri = "biojava";
    static final String CDATA = "CDATA";
    static final String PCDATA = "PCDATA";

    static String querySequenceType = null;
    static String hitSequenceType = null;

    /**
     * the SeqIOListener for this object
     */
    ContentHandler listener;

    /**
     *  Nesting class that provides callback interfaces to nested class
     */
    public StAXFeatureHandler staxenv;
    /**
     *  handler list for delegation
     */
    List handlers;

    /**
     *  are we in delegation from current object or is this pass my own?
     */
    int level = 0;

    /**
     *  This base class defines default behaviour for a StAX handler including
     *  delegation.
     */
    /**
     *  the version of constructor called depends on whether there is an
     *  explicit super(...) in the constructor of the derived class. If there
     *  is, that specific constructor will be called. If not, the parameterless
     *  will is called.
     */

    StAXFeatureHandler() {
        handlers = new ArrayList();
    }


    /**
     *  Constructor for the StAXFeatureHandler object
     *
     *@param  staxenv   Description of the Parameter
     */
    StAXFeatureHandler(StAXFeatureHandler staxenv) {
        handlers = new ArrayList();
        this.staxenv = staxenv;
    }

    // Class to implement bindings
    /**
     *  Description of the Class
     *
     *@author     david
     *@created    19 January 2002
     */
    class Binding {
        final ElementRecognizer recognizer;
        final StAXHandlerFactory handlerFactory;


        /**
         *  Constructor for the Binding object
         *
         *@param  er  Description of the Parameter
         *@param  hf  Description of the Parameter
         */
        Binding(ElementRecognizer er, StAXHandlerFactory hf) {
            recognizer = er;
            handlerFactory = hf;
        }
    }

    // method to add a handler
    // we do not distinguish whither it is a feature or property
    // handler.  The factory method creates the right type subclassed
    // from the correct type of handler
    /**
     *  Adds a feature to the Handler attribute of the StAXFeatureHandler object
     *
     *@param  rec      The feature to be added to the Handler attribute
     *@param  handler  The feature to be added to the Handler attribute
     */
    protected void addHandler(
            ElementRecognizer rec,
            StAXHandlerFactory handler) {
        handlers.add(new Binding(rec, handler));
    }

    /**
     * get the SeqIOListener for this parser
     */
    public ContentHandler getListener()
    {
//        System.out.println("in StAXFeatureHandler. staxenv is " + staxenv);
        return staxenv.listener;
    }

    /**
     *  Element-specific handler. Subclass this to do something useful!
     *
     *@param  nsURI             Description of the Parameter
     *@param  localName         Description of the Parameter
     *@param  qName             Description of the Parameter
     *@param  attrs             Description of the Parameter
     *@exception  SAXException  Description of the Exception
     */
    void startElementHandler(
            String nsURI,
            String localName,
            String qName,
            Attributes attrs)
             throws SAXException { }

    /**
     *  Handles basic entry processing for all feature handlers.
     *
     *@param  nsURI             Description of the Parameter
     *@param  localName         Description of the Parameter
     *@param  qName             Description of the Parameter
     *@param  attrs             Description of the Parameter
     *@param  dm                Description of the Parameter
     *@exception  SAXException  Description of the Exception
     */
    public void startElement(
            String nsURI,
            String localName,
            String qName,
            Attributes attrs,
            DelegationManager dm)
             throws SAXException {
        level++;

        // perform delegation
        // we must delegate only on features that are directly attached.
        // if I do not check that that's so, any element of a kind I delegate
        // on will be detected any depth within unrecognized tags.
        if (level == 2) {
//        System.out.println("StaxFeaturehandler.startElement starting. localName: " + localName + " " + level);
            for (int i = handlers.size() - 1; i >= 0; --i) {
                Binding b = (Binding) handlers.get(i);
                if (b.recognizer.filterStartElement(nsURI, localName, qName, attrs)) {
                    dm.delegate(b.handlerFactory.getHandler(staxenv));
                    return;
                }
            }
        }

        // call the element specific handler now.
        // remember that if we we have a delegation failure we pass here too!
        if (level == 1) {
            startElementHandler(nsURI, localName, qName, attrs);
        }
    }

    /**
     *  Element specific exit handler Subclass to do anything useful.
     *
     *@param  nsURI             Description of the Parameter
     *@param  localName         Description of the Parameter
     *@param  qName             Description of the Parameter
     *@param  handler           Description of the Parameter
     *@exception  SAXException  Description of the Exception
     */
    void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler handler)
             throws SAXException { }


    /**
     *  Handles basic exit processing.
     *
     *@param  nsURI             Description of the Parameter
     *@param  localName         Description of the Parameter
     *@param  qName             Description of the Parameter
     *@param  handler           Description of the Parameter
     *@exception  SAXException  Description of the Exception
     */
    public void endElement(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler handler)
             throws SAXException {
//      System.out.println("StAXFeatureHandler endElement called, localName, level: " + localName + " " + stackLevel);
        if ((--level) == 0) {
            endElementHandler(nsURI, localName, qName, handler);
        }
    }
}

