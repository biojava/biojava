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
import java.util.StringTokenizer;
import java.util.Vector;

import org.biojava.bio.BioException;
import org.biojava.utils.ParseErrorEvent;
import org.biojava.utils.ParseErrorListener;
import org.biojava.utils.ParseErrorSource;

/**
 * Simple filter which handles attribute lines from an EMBL file. This
 * class delegates creation of <code>Feature</code>s to a
 * <code>FeatureTableParser</code>, which in turn delegates creation
 * of <code>Locations</code> to an <code>EmblLikeLocationParser</code>
 * which is shared with the <code>GenbankProcessor</code>.
 *
 * An <code>EmblLikeLocationParser</code> parses EMBL/Genbank style
 * locations. Supported location forms:
 *
 * <pre>
 *   123
 *  <123 or >123
 *  (123.567)
 *  (123.567)..789
 *   123..(567.789)
 *  (123.345)..(567.789)
 *   123..456
 *  <123..567 or 123..>567 or <123..>567
 *   123^567
 *   AL123465:(123..567)
 * </pre>
 *
 * The only EMBL header information retained over a read/write cycle
 * is the accession number (all numbers).
 *
 * @author Thomas Down
 * @author Greg Cox
 * @author Keith James
 * @since 1.1
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */

public class EmblProcessor
	extends
		SequenceBuilderFilter
	implements
		ParseErrorSource
{
    public static final String PROPERTY_EMBL_ACCESSIONS = "embl_accessions";

    private boolean mBadFeature = false;
    private Vector mListeners = new Vector();

    /**
     * Factory which wraps SequenceBuilders in an EmblProcessor
     *
     * @author Thomas Down
     */

    public static class Factory implements SequenceBuilderFactory, Serializable {
	private SequenceBuilderFactory delegateFactory;

	public Factory(SequenceBuilderFactory delegateFactory) {
	    this.delegateFactory = delegateFactory;
	}

	public SequenceBuilder makeSequenceBuilder() {
	    return new EmblProcessor(delegateFactory.makeSequenceBuilder());
	}
    }

    private FeatureTableParser features;

    public EmblProcessor(SequenceBuilder delegate) {
        super(delegate);
        features = new FeatureTableParser(this, "EMBL");
    }

    public void endSequence() throws ParseException {
        // Avoids leaving a null name and null URI if there is no
        // accession number. If accession number is vital, failure of
        // test of accessions.size() > 0 should throw a
        // ParseException.
        //String  id = "";
        String uri = "";
        if (accessions.size() > 0) {
            //id = (String) accessions.get(0);
            uri = "urn:sequence/embl:" + (String) accessions.get(0);
            getDelegate().addSequenceProperty(PROPERTY_EMBL_ACCESSIONS, accessions);
        }

        //getDelegate().setName(id);
        getDelegate().setURI(uri);
        getDelegate().endSequence();
    }

    private List accessions;

    {
        accessions = new ArrayList();
    }

    public void addSequenceProperty(Object key, Object value)
        throws ParseException
    {
        try
        {
            if (mBadFeature)
            {
                // If this feature is bad in some way, ignore it.
                if (value != null)
                {
                    String featureLine = value.toString();
                    if((key.equals(EmblLikeFormat.FEATURE_TABLE_TAG)) &&
                       (featureLine.charAt(0) != ' '))
                    {
                        // If the offending feature is past, start reading data again
                        mBadFeature = false;
                        features.startFeature(featureLine.substring(0, 15).trim());
                        features.featureData(featureLine.substring(16));
                    }
                }
            }
            else
            {
                // Tidy up any end-of-block jobbies
                if (features.inFeature() &&
                    !key.equals(EmblLikeFormat.FEATURE_TABLE_TAG))
                {
                    features.endFeature();
                }

                if (key.equals(EmblLikeFormat.FEATURE_TABLE_TAG))
                {
                    String featureLine = value.toString();
                    if (featureLine.charAt(0) != ' ')
                    {
                        // This is a featuretype field
                        if (features.inFeature())
                        {
                            features.endFeature();
                        }

                        features.startFeature(featureLine.substring(0, 15).trim());
                    }
                    features.featureData(featureLine.substring(16));
                }
                else
                {
                    getDelegate().addSequenceProperty(key, value);

                    if (key.equals(EmblLikeFormat.ACCESSION_TAG))
                    {
                        String acc = value.toString();
                        StringTokenizer toke = new StringTokenizer(acc, "; ");
                        while (toke.hasMoreTokens())
                        {
                            accessions.add(toke.nextToken());
                        }
                    }
                    else if (key.equals(EmblLikeFormat.ID_TAG)) {
                        StringTokenizer toke = new StringTokenizer((String) value);
                        getDelegate().setName(toke.nextToken());
                    }
                }
            }
        }
        catch (BioException ex)
        {
            // If an exception is thrown, read past the offending feature
            mBadFeature = true;
            ParseErrorEvent offendingLineEvent =
                new ParseErrorEvent(this, "This line could not be parsed: "
                                    + value.toString());
            this.notifyParseErrorEvent(offendingLineEvent);
        }
        catch (IndexOutOfBoundsException ex)
        {
            // This occurs when for some line min > max
            mBadFeature = true;
            ParseErrorEvent offendingLineEvent =
                new ParseErrorEvent(this, "From must be less than To: "
                                    + value.toString());
            this.notifyParseErrorEvent(offendingLineEvent);
        }
    }

    /**
     * Adds a parse error listener to the list of listeners if it isn't already
     * included.
     *
     * @param theListener Listener to be added.
     */
    public synchronized void addParseErrorListener(ParseErrorListener theListener)
    {
        if (mListeners.contains(theListener) == false)
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
    public synchronized void removeParseErrorListener(ParseErrorListener theListener)
    {
        if (mListeners.contains(theListener) == true)
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
