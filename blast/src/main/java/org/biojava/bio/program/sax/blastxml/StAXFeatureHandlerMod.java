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

import org.biojava.utils.stax.DelegationManager;
import org.xml.sax.Attributes;
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

class StAXFeatureHandlerMod extends StAXFeatureHandler {

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
        if (level == 1) {
//        System.out.println("StaxFeaturehandler.startElement starting. localName: " + localName + " " + level);
            for (int i = handlers.size() - 1; i >= 0; --i) {
                Binding b = (Binding) handlers.get(i);
                if (b.recognizer.filterStartElement(nsURI, localName, qName, attrs)) {
                    dm.delegate(b.handlerFactory.getHandler(this));
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
}

