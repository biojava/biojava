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

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * StAX handler for GAME &lt;residues&gt; elements.
 * derived from Thomas Down's PropDetailHandler
 *
 * <p>
 * This takes the sequence supplied by &lt;residues&gt; elements
 * and feeds it to a StreamParser associated with a SeqIOLIstener
 * instance.
 *
 * @author David Huen
 * @author Thomas Down
 * @since 1.8
 */
public class GAMEResiduesPropHandler extends SequenceContentHandlerBase {
  public static final StAXHandlerFactory GAME_RESIDUES_PROP_HANDLER_FACTORY = new StAXHandlerFactory() {
	  public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                 return new GAMEResiduesPropHandler(staxenv);
	  }
            } ;

  private StAXFeatureHandler staxenv;

  public GAMEResiduesPropHandler(StAXFeatureHandler staxenv) {
    super();
    this.staxenv = staxenv;
  }

  public void startElement(String nsURI,
                                    String localName,
                                    String qName,
                                    Attributes attrs,
                                    DelegationManager dm)
                     throws SAXException
  {
    super.startElement(nsURI, localName, qName, attrs, dm);

    // set up StreamParser
    try {
	SymbolTokenization tokens = DNATools.getDNA().getTokenization("token");
	super.setStreamParser(tokens.parseStream(staxenv.featureListener));
    } catch (BioException ex) {
	throw new BioError("Assertion failed: couldn't get tokenization from DNA alphabet");
    }
  }

}
