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
import java.util.StringTokenizer;

/**
 * Simple filter which performs a default extraction of data from
 * the description lines of FASTA files.  Behaviour is similar
 * to DefaultDescriptionReader in the old I/O framework.
 *
 * @author Thomas Down
 * @since 1.1
 * @deprecated Use org.biojavax.bio.seq.io.FastaFormat
 */

public class FastaDescriptionLineParser extends SequenceBuilderFilter {
    /**
     * Factory which wraps SequenceBuilders in a FastaDescriptionLineParser
     *
     * @author Thomas Down
     */

    public static class Factory implements SequenceBuilderFactory, Serializable {
	private SequenceBuilderFactory delegateFactory;

	public Factory(SequenceBuilderFactory delegateFactory) {
	    this.delegateFactory = delegateFactory;
	}

	public SequenceBuilder makeSequenceBuilder() {
	    return new FastaDescriptionLineParser(delegateFactory.makeSequenceBuilder());
	}
    }

    public FastaDescriptionLineParser(SequenceBuilder delegate) {
	super(delegate);
    }

    public void addSequenceProperty(Object key, Object value) throws ParseException {
	getDelegate().addSequenceProperty(key, value);

	if (FastaFormat.PROPERTY_DESCRIPTIONLINE.equals(key)) {
	    String dline = value.toString();
	    StringTokenizer toke = new StringTokenizer(dline);
	    String name = toke.nextToken();
	    setName(name);
	    setURI("urn:sequence/fasta:" + name);
	    if (toke.hasMoreTokens()) {
		getDelegate().addSequenceProperty("description", toke.nextToken("******"));
	    }
	} 
    }
}
