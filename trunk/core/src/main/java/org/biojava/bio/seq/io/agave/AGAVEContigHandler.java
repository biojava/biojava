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
import java.util.List;
import java.util.ListIterator;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleAssembly;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *
 * Handles the AGAVE &lt;contig&gt; element
 * @author Hanning Ni    Doubletwist Inc
 * @author Greg Cox
 */
public class AGAVEContigHandler
               extends StAXFeatureHandler  implements AGAVEContigCallbackItf, SequenceHandler
{
  public static final StAXHandlerFactory AGAVE_CONTIG_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEContigHandler(staxenv);
    }
  };

  protected Sequence sequence ;
  private List sequenceSet ;
  private SymbolList dna ;

  AGAVEContigHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("contig", true);
    sequenceSet = new ArrayList(1) ;

    // setup handlers
       // <db_id>
       super.addHandler(new ElementRecognizer.ByLocalName("db_id"),
         AGAVEDbIdPropHandler.AGAVE_DBID_PROP_HANDLER_FACTORY);
      //view
       super.addHandler(new ElementRecognizer.ByLocalName("view"),
         AGAVEViewPropHandler.AGAVE_VIEW_PROP_HANDLER_FACTORY);
       // <note>
       super.addHandler(new ElementRecognizer.ByLocalName("note"),
         AGAVENotePropHandler.AGAVE_NOTE_PROP_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("fragment_order"),
         AGAVEFragmentOrderHandler.AGAVE_FRAGMENT_ORDER_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("unordered_fragments"),
         AGAVEUnorderedFragmentsHandler.AGAVE_UNORDERED_FRAGMENTS_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("assembly"),
         AGAVEAssemblyHandler.AGAVE_ASSEMBLY_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("sequence"),
         AGAVESeqPropHandler.AGAVE_SEQ_PROP_HANDLER_FACTORY);
       //<sequence_map>
       super.addHandler(new ElementRecognizer.ByLocalName("sequence_map"),
         AGAVESeqMapHandler.AGAVE_SEQ_MAP_HANDLER_FACTORY);
       //<map_location>
       super.addHandler(new ElementRecognizer.ByLocalName("map_location"),
         AGAVEMapLocationPropHandler.AGAVE_MAP_LOCATION_PROP_HANDLER_FACTORY);
  }
  public void reportSequence(Sequence sequence)
  {
      sequenceSet.add( sequence ) ;
  }

  public void reportDna(String dna)
  {
       try{
       StringBuffer sb = new StringBuffer()  ;
       for( int i = 0 ; i < dna.length(); i++)
       {
           char c = dna.charAt(i) ;
           if( c != ' '  && c != '\n' && c!= '\t')
              sb.append( c  );
       }
     this.dna = DNATools.createDNA( sb.substring(0) );
     }catch(Exception e){
       e.printStackTrace( );
     }
  }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
                throws SAXException
  {
     try{
      featureListener.startSequence();
      boolean forFeature = true ;
      setProperty( "length",  attrs.getValue("length") , forFeature) ;
      }catch(Exception e){
         throw new SAXException( e.getMessage() ) ;
      }
  }

   /*
   protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template st = new StrandedFeature.Template();

    // assume feature set to describe a transcript
    st.type = "contig";
    st.strand = StrandedFeature.UNKNOWN;
    // set up annotation bundle
    st.annotation = annot; //new SmallAnnotation();
    st.location = new  Location.EmptyLocation();
    if( staxenv != null )
        staxenv. subFeatures .add( this ) ;

    return st;
  }*/

  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
                throws SAXException
  {
       try{
       if( sequenceSet.size() == 0 )
       {
            if( dna  == null )
                throw new SAXException("dna sequence must offered ") ;
            sequence  = new SimpleSequence( dna, " ", "simple_sequence " , annot ) ;
            //Feature feature = sequence.createFeature( featureTemplate ) ;
            //realizeSubFeatures( feature ) ;
            addFeatureToSequence(sequence) ;
       }
       else
       {
           int contig_len =new Integer( (String) UtilHelper.getProperty(featureTemplate.annotation, "length" )).intValue()  ;
           sequence = new SimpleAssembly(contig_len , "contig", "contig") ;
         //  sequence.createFeature( featureTemplate ) ;
           ComponentFeature.Template cft = new ComponentFeature.Template();
           int global_start = 1 ;
           for (ListIterator i = sequenceSet.listIterator(); i.hasNext();)
           {
               Sequence seq = ( Sequence ) i.next() ;
             //   ComponentFeature.Template cft = new ComponentFeature.Template();
                cft.type = "fragment";
                cft.source = "contig";
                cft.annotation = Annotation.EMPTY_ANNOTATION;
                cft.strand = StrandedFeature.POSITIVE;
                cft.location = new RangeLocation(global_start, global_start +  seq.length() - 1 );
                cft.componentSequence = seq;
                cft.componentLocation = new RangeLocation(1, seq.length());
                sequence.createFeature(cft);
                global_start += seq.length() ;
           }
       }
       //add sequence to chromsome
       appendToTop(sequence, staxenv) ;
       featureListener.endSequence() ;
       }catch(Exception e){
         throw new SAXException( e.getMessage() ) ;
       }
    }
    private void appendToTop(Sequence sequence, StAXFeatureHandler staxenv)
    {
              if (staxenv instanceof AGAVEChromosomeCallbackItf)
                {
                    ((AGAVEChromosomeCallbackItf) staxenv).reportSequence( sequence );
                    return;
                }
                 if (staxenv instanceof AGAVECallbackItf)
                {
                    ((AGAVECallbackItf) staxenv).reportSequence( sequence );
                    return;
                }
               appendToTop(sequence, staxenv.staxenv );
    }

}

