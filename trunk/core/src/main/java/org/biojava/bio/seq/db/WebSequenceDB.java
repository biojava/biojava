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

package org.biojava.bio.seq.db;

import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.seq.io.StreamReader;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeVetoException;

/**
 * Functions for access to a web based database that returns sequences
 * in a variety of formats.
 *
 * @author Jason Stajich
 * @author Matthew Pocock
 * @author Mark Schreiber
 * @author Richard Holland
 */

public abstract class WebSequenceDB
extends AbstractChangeable
implements SequenceDBLite {
  protected abstract SequenceFormat getSequenceFormat();

  protected abstract URL getAddress(String id)
  throws MalformedURLException;

  protected abstract Alphabet getAlphabet();

  /**
   * Gets a sequence using its unique ID (eg for GenBank this would be the GI number)
   * @param id the unique ID
   * @return the matching sequence
   * @throws BioException if the ID is invalid
   * @throws BioException if the io operation times out or has problems
   *    connecting. Can also indicate an invalid URL has been constructed.
   */
  public Sequence getSequence(String id)
  throws BioException {
    if( id.equals("") ) {
      throw new BioException("did not specify a valid id for getSequence");
    }

    try {
      URL queryURL = getAddress(id);
      //System.err.println("query is "+ queryURL.toString());
      URLConnection connection = queryURL.openConnection();
      SequenceFormat sFormat = getSequenceFormat();

//      SequenceBuilder sbuilder = new SimpleSequenceBuilder();
//      FastaDescriptionLineParser sFact =
//        new FastaDescriptionLineParser(sbuilder);

      Alphabet alpha = getAlphabet();
      SequenceBuilderFactory sFact = SeqIOTools.formatToFactory(sFormat,alpha);
      SymbolTokenization rParser = alpha.getTokenization("token");
      //System.err.println("got data from "+ queryURL);
      SequenceIterator seqI = new StreamReader(
        connection.getInputStream(),
        sFormat, rParser, sFact
      );

      return seqI.nextSequence();
    } catch ( Exception e ){
      throw new BioException(e);
    }
  }

  /**
   * Not supported, You can't add sequences to a WebDB!
   * @param seq the sequence you tried to add
   * @throws ChangeVetoException always!
   */
  public void addSequence(Sequence seq)
  throws ChangeVetoException {
    throw new ChangeVetoException(
      "Can't add sequences from web sequence DB: " +
      seq.getName()
    );
  }

  /**
   * Not supported, you can't remove a sequence from a WebDB!
   * @param id the sequence you tried to change.
   * @throws ChangeVetoException always!
   */
  public void removeSequence(String id)
  throws ChangeVetoException {
    throw new ChangeVetoException(
      "Can't remove sequences from web sequence DB: " +
      id
    );
  }
}
