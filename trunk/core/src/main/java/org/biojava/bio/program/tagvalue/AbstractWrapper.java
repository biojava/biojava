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


package org.biojava.bio.program.tagvalue;

import org.biojava.utils.ParserException;

/**
 * <p>
 * An abstract TagValueWrapper that does nothing!
  </p>
 * <p>
 * Useful for writing TagValueWrappers which
 * only act on a subset of the events.
 * </p>
 * @author David Huen
 */
public class AbstractWrapper
	implements TagValueWrapper
{
	private TagValueListener delegate;

	public void setDelegate(TagValueListener delegate)
	{
		this.delegate = delegate;
	}

	public TagValueListener getDelegate()
	{
		return delegate;
	}

	public void startRecord() throws ParserException {}
	public void startTag(Object tag) throws ParserException {}
	public void endTag() throws ParserException {}
	public void endRecord() throws ParserException {}
	public void value(TagValueContext ctxt, java.lang.Object value) throws ParserException {}
}
