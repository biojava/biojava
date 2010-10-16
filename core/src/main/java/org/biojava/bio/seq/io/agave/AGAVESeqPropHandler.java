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

package org.biojava.bio.seq.io.agave;

import java.util.ListIterator;

import org.xml.sax.SAXException;





/**

 * Deals with sequence code

 * 

 * @author Hanning Ni    Doubletwist Inc
  * @author Greg Cox

 */

public class AGAVESeqPropHandler

               extends StAXPropertyHandler

{

  // set up factory method

  public static final StAXHandlerFactory AGAVE_SEQ_PROP_HANDLER_FACTORY

    = new StAXHandlerFactory() {

    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {

      return new AGAVESeqPropHandler(staxenv);

    }

  };

  private StringBuffer dnaTokens  ;



  AGAVESeqPropHandler(StAXFeatureHandler staxenv) {

    // execute superclass method to setup environment

    super(staxenv);

    setHandlerCharacteristics("sequence", true);

  }



   public void characters(char[] ch, int start, int length)

        throws SAXException

  {

       dnaTokens = new StringBuffer()  ;

       for( int i = start ; i < start + length; i++)

       {

           char c = ch[i] ;

           if( c != ' '  && c != '\n' && c!= '\t')

              dnaTokens.append( c  );

       }

  }



  public void endElementHandler(

                String nsURI,

                String localName,

                String qName,

                StAXContentHandler handler)

              throws SAXException

  {

        int currLevel = staxenv.getLevel();

        if (currLevel >=1)

        {

            ListIterator li = staxenv.getHandlerStackIterator(currLevel);

            while(li.hasPrevious())

            {

               Object ob = li.previous();

                if (ob instanceof AGAVEBioSeqCallbackItf)

                {

                     ((AGAVEBioSeqCallbackItf) ob).reportDna( dnaTokens.substring(0) );

                    return;

                }

            }

            li = staxenv.getHandlerStackIterator(currLevel);

            while (li.hasPrevious())

            {

                Object ob = li.previous();

                if (ob instanceof AGAVEContigCallbackItf)

                {

                    ((AGAVEContigCallbackItf) ob).reportDna( dnaTokens.substring(0) );

                    return;

                }

            }

        }

  }

}

