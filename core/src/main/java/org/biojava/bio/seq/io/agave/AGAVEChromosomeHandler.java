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
package org.biojava.bio.seq.io.agave;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import org.biojava.bio.seq.Sequence;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *
 * Handles the AGAVE &lt;chromosome&gt; element
 *
 * @author Hanning Ni     Doubletwist Inc
 */
public class AGAVEChromosomeHandler
               extends StAXFeatureHandler  implements AGAVEChromosomeCallbackItf, SequenceHandler
{
  public static final StAXHandlerFactory AGAVE_CHROMOSOME_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEChromosomeHandler();
    }
  };

  private List sequenceSet ;
  private String chromNum ;
  private String chromLen ;

  AGAVEChromosomeHandler() {
    super();
    setHandlerCharacteristics("chromosome", true);
    sequenceSet = new ArrayList(1) ;

    // setup handlers
       // ignore it , view
       super.addHandler(new ElementRecognizer.ByLocalName("view"),
         AGAVEViewPropHandler.AGAVE_VIEW_PROP_HANDLER_FACTORY);
      // <note>
       super.addHandler(new ElementRecognizer.ByLocalName("contig"),
         AGAVEContigHandler.AGAVE_CONTIG_HANDLER_FACTORY);
  }

  public void reportSequence(Sequence sequence)
  {
      sequenceSet.add( sequence ) ;
  }
  public Iterator getSequences()
  {
      return sequenceSet.iterator() ;
  }
  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
                throws SAXException
  {
       chromNum = attrs.getValue("number") ;
       chromLen = attrs.getValue("length") ;
  }

  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
                throws SAXException
  {
       try{
       for (ListIterator i = sequenceSet.listIterator(); i.hasNext();)
       {
           Sequence seq = ( Sequence ) i.next() ;
           if( chromNum !=  null )
               seq. getAnnotation().setProperty( "chromosome_number", chromNum);
           if( chromLen != null )
               seq. getAnnotation().setProperty( "chromosome_length", chromLen);

           appendToTop(seq, staxenv) ;
       }
       }catch(Exception e){
          throw new SAXException( e.getMessage() ) ;
        }

  }

     private void appendToTop(Sequence sequence, StAXFeatureHandler staxenv)
    {
        if( staxenv instanceof AGAVECallbackItf)
        {
             ((AGAVECallbackItf) staxenv).reportSequence( sequence );
             return;
        }
        else appendToTop(sequence, staxenv.staxenv );
    }

}

