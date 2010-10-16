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

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.impl.DummySequence;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeVetoException;

/**
 * <code>DummySequenceDB</code> is an implementation which contains
 * only a <code>DummySequence</code>. It will return the same
 * <code>DummySequence</code> instance regardless of the sequence id
 * used to retrieve a sequence.
 *
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a>
 * @since 1.2
 */
public class DummySequenceDB extends AbstractSequenceDB
    implements SequenceDB
{
    private String name;
    private Map    seqs;

    public DummySequenceDB(String name)
    {
        this.name = name;
        seqs = new HashMap();
        seqs.put("dummy", new DummySequence("dummy", "dummy"));

        // Lock to prevent sequences being added or removed
        this.addChangeListener(ChangeListener.ALWAYS_VETO);
    }

    public Set ids()
    {
        return Collections.unmodifiableSet(seqs.keySet());
    }

    public SequenceIterator sequenceIterator()
    {
        return new SequenceIterator()
            {
                Iterator i = seqs.values().iterator();

                public boolean hasNext()
                {
                    return i.hasNext();
                }

                public Sequence nextSequence()
                {
                    return (Sequence) i.next();
                }
            };
    }

    public String getName()
    {
        return name;
    }

    public Sequence getSequence(String id)
        throws IllegalIDException, BioException
    {
        return (Sequence) seqs.get("dummy");
    }

    public void addSequence(Sequence seq)
        throws IllegalIDException, BioException, ChangeVetoException { }

    public void removeSequence(String id)
        throws IllegalIDException, BioException, ChangeVetoException { }
}
