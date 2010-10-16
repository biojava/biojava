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

package org.biojava.bio.seq.io;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import org.biojava.bio.BioException;
import org.biojava.utils.ParseErrorEvent;
import org.biojava.utils.ParseErrorListener;
import org.biojava.utils.ParseErrorSource;

/**
 * Simple filter which handles attribute lines from a Genbank file
 *
 * @author Greg Cox
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */

public class GenbankProcessor extends SequenceBuilderFilter
implements ParseErrorSource
{
  public static final String PROPERTY_GENBANK_ACCESSIONS = "genbank_accessions";
  private boolean mBadFeature = false;
  private Vector mListeners = new Vector();
  
  /**
  * Factory which wraps sequence builders in a GenbankProcessor
  *
  * @author Greg Cox
  */
  public static class Factory implements SequenceBuilderFactory, Serializable
  {
    private SequenceBuilderFactory delegateFactory;
    
    public Factory(SequenceBuilderFactory theDelegate)
    {
      delegateFactory = theDelegate;
    }
    
    public SequenceBuilder makeSequenceBuilder()
    {
      return new GenbankProcessor(delegateFactory.makeSequenceBuilder());
    }
    
    public SequenceBuilder makeSequenceBuilder(String theSource)
    {
      return new GenbankProcessor(delegateFactory.makeSequenceBuilder(), theSource);
    }
  }
  
  protected FeatureTableParser features;
  private List accessions;

  {
    accessions = new ArrayList();
  }
  
  public GenbankProcessor(SequenceBuilder theDelegate, String theSource)
  {
    super(theDelegate);
    features = new FeatureTableParser(this, theSource);
  }
  
  public GenbankProcessor(SequenceBuilder theDelegate)
  {
    super(theDelegate);
    features = new FeatureTableParser(this, "GenBank");
  }
  
  public void endSequence() throws ParseException
  {
    // Avoids leaving a null name and null URI if there is no
    // accession number. If accession number is vital, failure of
    // test of accessions.size() > 0 should throw a
    // ParseException.
    String uri = "";
    if (accessions.size() > 0) {
      uri = "urn:sequence/genbank:" + (String) accessions.get(0);
      getDelegate().addSequenceProperty(PROPERTY_GENBANK_ACCESSIONS, accessions);
    }
    
    getDelegate().setURI(uri);
    getDelegate().endSequence();
  }
  
  public void addSequenceProperty(Object key, Object value) throws ParseException
  {
    try
    {
      if(mBadFeature)
      {
        // If this feature is bad in some way, ignore it.
        String featureLine = value.toString();
        if((key.equals(GenbankFormat.FEATURE_FLAG)) && (featureLine.charAt(0) != ' '))
        {
          // If the offending feature is past, start reading data again
          mBadFeature = false;
          features.startFeature(featureLine.substring(0, 16).trim());
          features.featureData(featureLine.substring(16));
        }
      }
      else
      {
        if (features.inFeature() && !(key.equals(GenbankFormat.FEATURE_FLAG)))
        {
          features.endFeature();
        }
        
        if(key.equals(GenbankFormat.FEATURE_FLAG))
        {
          String featureLine = value.toString();
          if (featureLine.charAt(0) != ' ')
          {
            // This is a featuretype field
            if (features.inFeature())
            {
              features.endFeature();
            }
            features.startFeature(featureLine.substring(0, 16).trim());
          }
          features.featureData(featureLine.substring(16));
        }
        else
        {
          getDelegate().addSequenceProperty(key, value);
          if (key.equals(GenbankFormat.ACCESSION_TAG))
          {
            accessions.add(value);
          }
          else if (key.equals(GenbankFormat.LOCUS_TAG)) {
            getDelegate().setName((String) value);
            features.setSeqID((String) value);
          }
        }
      }
    }
    catch (BioException ex)
    {
      // If an exception is thrown, notify the listeners and read past
      // the offending feature
      mBadFeature = true;
      ParseErrorEvent offendingLineEvent = new ParseErrorEvent(this, "This line could not be parsed: " + value.toString());
      this.notifyParseErrorEvent(offendingLineEvent);
    }
    catch (IndexOutOfBoundsException ex)
    {
      // This occurs when for some line min > max
      mBadFeature = true;
      ParseErrorEvent offendingLineEvent = new ParseErrorEvent(this, "From must be less than To: " + value.toString());
      this.notifyParseErrorEvent(offendingLineEvent);
    }
  }
  
  /**
  * Adds a parse error listener to the list of listeners if it isn't already
  * included.
  *
  * @param theListener Listener to be added.
  */
  public synchronized void addParseErrorListener(
  ParseErrorListener theListener)
  {
    if(mListeners.contains(theListener) == false)
    {
      mListeners.addElement(theListener);
    }
  }
  
  /**
  * Removes a parse error listener from the list of listeners if it is
  * included.
  *
  * @param theListener Listener to be removed.
  */
  public synchronized void removeParseErrorListener(
  ParseErrorListener theListener)
  {
    if(mListeners.contains(theListener) == true)
    {
      mListeners.removeElement(theListener);
    }
  }
  
  // Protected methods
  /**
  * Passes the event on to all the listeners registered for ParseErrorEvents.
  *
  * @param theEvent The event to be handed to the listeners.
  */
  protected void notifyParseErrorEvent(ParseErrorEvent theEvent)
  {
    Vector listeners;
    synchronized(this)
    {
      listeners = (Vector)mListeners.clone();
    }
    
    for (int index = 0; index < listeners.size(); index++)
    {
      ParseErrorListener client = (ParseErrorListener)listeners.elementAt(index);
      client.BadLineParsed(theEvent);
    }
  }
}
