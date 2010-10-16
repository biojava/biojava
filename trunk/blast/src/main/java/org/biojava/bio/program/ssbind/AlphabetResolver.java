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

package org.biojava.bio.program.ssbind;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.FiniteAlphabet;

/**
 * <code>AlphabetResolver</code>s are helpers which determine which
 * type of sequence <code>Alphabet</code> to expect from a search
 * result. Now public to allow use by anyone making custom handlers.
 *
 * @author Keith James
 * @since 1.2
 */
public class AlphabetResolver
{
    static final int     DNA = 0;
    static final int PROTEIN = 1;

    /**
     * <code>resolveAlphabet</code> returns an appropriate
     * <code>Alphabet</code> for an arbitrary identifier. The protein
     * alphabet returned will include the termination character as
     * e.g. BLASTX 6-frame translations are likely to include stops.
     *
     * @param identifier a <code>String</code> identifier (recognised
     * are BLASTN, BLASTP, BLASTX, TBLASTN, TBLASTX, DNA and PROTEIN).
     *
     * @return a <code>FiniteAlphabet</code>.
     *
     * @exception BioException if the identifier is not known.
     */
    public static FiniteAlphabet resolveAlphabet(String identifier)
        throws BioException
    {
        int type = 0;

        identifier = identifier.toUpperCase();

        if (identifier.indexOf("TBLASTN") != -1)
            type = PROTEIN;
        else if (identifier.indexOf("TBLASTX") != -1)
            type = PROTEIN;
        else if (identifier.indexOf("BLASTN") != -1)
            type = DNA;
        else if (identifier.indexOf("BLASTP") != -1)
            type = PROTEIN;
        else if (identifier.indexOf("BLASTX") != -1)
            type = PROTEIN;
        else if (identifier.indexOf("DNA") != -1)
            type = DNA;
        else if (identifier.indexOf("PROTEIN") != -1)
            type = PROTEIN;
        else
            throw new BioException("Failed to resolve sequence type from identifier '"
                                   + identifier
                                   + "'");

        switch (type)
        {
            case DNA:
                return DNATools.getDNA();

            case PROTEIN:
                return ProteinTools.getTAlphabet();

            default:
                throw new BioError("Internal error in AlphabetResolver: failed to resolve to either DNA or protein alphabets");
        }
    }
}
