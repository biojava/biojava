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

import org.biojava.bio.seq.Sequence;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * Handles the root AGAVE element
 * modified for agave format
 * 
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVEHandler extends StAXFeatureHandler implements AGAVECallbackItf{

  // there is only one AGAVE element encompassing the entire file
  // AGAVE files.

  //store the sequences from each direct-subtag of sciobj
  private List sequenceSet ;

  public AGAVEHandler() {
     super() ;
     sequenceSet = new ArrayList(1) ;
     setHandlerCharacteristics("sciobj", true);

    // setup handlers
       // <bio_sequence>
       super.addHandler(new ElementRecognizer.ByLocalName("bio_sequence"),
         AGAVEBioSeqHandler.AGAVE_BIO_SEQ_HANDLER_FACTORY);
       // <contig>
       super.addHandler(new ElementRecognizer.ByLocalName("contig"),
         AGAVEContigHandler.AGAVE_CONTIG_HANDLER_FACTORY);
       // <computation> not used and handled yet
        super.addHandler(new ElementRecognizer.ByLocalName("computation"),
        AGAVEComputationHandler.AGAVE_COMPUTATION_HANDLER_FACTORY);
       //chromosome
       super.addHandler(new ElementRecognizer.ByLocalName("chromosome"),
         AGAVEChromosomeHandler.AGAVE_CHROMOSOME_HANDLER_FACTORY);
  }

  /**
   * @param sequence from sub-tag &lt;bio_sequence&gt;/&lt;contig&gt;/&lt;chromosome&gt;
   * <pre>
   * bio_sequence --> SimpleSequence
   * contig   --> SimpleAssembly
   *          --> SimpleSequence( if only one sequence )
   * chromosome -> SimpleAssembly
   *            -> SimpleSequence( if only one sequence)
   * </pre>
   */
  public void reportSequence(Sequence sequence)
  {
      sequenceSet.add( sequence ) ;
  }

  /**
   * get all the top level sequences
   * bio_sequence --> SimpleSequence
   * contig   --> SimpleAssembly
   *          --> SimpleSequence( if only one sequence )
   * chromosome -> SimpleAssembly
   *            -> SimpleSequence( if only one sequence)
   */
  public Iterator getSequences()
  {
      return sequenceSet.iterator() ;
  }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs,
                DelegationManager dm)
         throws SAXException
  {
    // check the element type
    if (!(localName.equals( "sciobj")) )
      throw new SAXException("root element of file is not a sciobj element");

    // check file version
    String version = attrs.getValue("version");
    if (! version.startsWith("2")  )
      throw new SAXException("AGAVE version is not 2.*,  we only support 2.* ");
 }


}

