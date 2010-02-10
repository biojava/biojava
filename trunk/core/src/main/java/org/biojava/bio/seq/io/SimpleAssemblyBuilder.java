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

package org.biojava.bio.seq.io;

import java.io.NotSerializableException;
import java.io.ObjectStreamException;
import java.io.Serializable;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleAssembly;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.StaticMemberPlaceHolder;

/**
 * Basic SequenceBuilder implementation which accumulates all
 * notified information and creates a SimpleAssembly.
 *
 * @author David Huen
 * @version 1.2
 */

public class SimpleAssemblyBuilder extends SequenceBuilderBase {
    public final static SequenceBuilderFactory FACTORY = new SSBFactory();

    private static class SSBFactory implements SequenceBuilderFactory, Serializable {
	private SSBFactory() {
	}

	public SequenceBuilder makeSequenceBuilder() {
	    return new SimpleAssemblyBuilder();
	}

	private Object writeReplace() throws ObjectStreamException {
	    try {
		return new StaticMemberPlaceHolder(SimpleAssemblyBuilder.class.getField("FACTORY"));
	    } catch (NoSuchFieldException nsfe) {
		throw new NotSerializableException(nsfe.getMessage());
	    }
	}
    }

    private void checkSeq()
    {
      // check that seq exists: if not, create it.
      // this is done to permit lazy instatiation of the SimpleAssembly
      // which gives time for name and uri to be set.
      if (seq == null) seq = new SimpleAssembly(name, uri);
    }

    //
    // SeqIOListener
    //

    public void addSymbols(Alphabet alpha, Symbol[] syms, int pos, int len)
        throws IllegalAlphabetException
    {
      System.err.println("SimpleAssemblyBuilder: illegal attempt to add symbols");
    }

    public ComponentFeature addComponentSequence(ComponentFeature.Template cft)
      throws BioException, ChangeVetoException
    {
        checkSeq();

        return (ComponentFeature) seq.createFeature(cft);
    }

  public Sequence makeSequence()
          throws BioException
  {
    checkSeq();

    return super.makeSequence();
  }
}
