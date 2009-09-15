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

/**
 * Class description
 *
 * @author Greg Cox
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */
public class ProteinRefSeqProcessor extends GenbankProcessor
{
// Static variables

// Member variables

// Constructors and initialization
	public ProteinRefSeqProcessor(SequenceBuilder theDelegate)
	{
		super(theDelegate);
		features = new FeatureTableParser(this, "RefSeq:Protein");
	}

// Interface implementations

// Public methods
	/**
	 * Factory which wraps sequence builders in a ProteinRefSeqProcessor
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
			return new ProteinRefSeqProcessor(delegateFactory.makeSequenceBuilder());
		}
	}

// Protected methods

// Private methods
}
