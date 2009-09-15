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

/**
 * <code>AlignIOConstants</code> contains constants used to identify
 * sequence formats, alphabets etc, in the context of reading and
 * writing alignments.
 *
 * <p>An <code>int</code> used to specify symbol alphabet and
 * sequence format type is derived thus:</p>
 *
 * <ul>
 *   <li>
 *    The two least significant bytes are reserved for format types
 *    such as MSF, CLUSTAL etc.
 *   </li>
 *
 *   <li>
 *    The two most significant bytes are reserved for alphabet and
 *    symbol information such as AMBIGUOUS, DNA, RNA, AA etc.
 *   </li>
 *
 *   <li>
 *    Bitwise OR combinations of each component <code>int</code> are used
 *    to specify combinations of format type and symbol information. To
 *    derive an <code>int</code> identifier for DNA with ambiguity codes
 *    in Fasta format, bitwise OR the AMBIGUOUS, DNA and FASTA values.
 *   </li>
 * </ul>
 *
 * @author Keith James
 */
public final class AlignIOConstants
{
    /**
     * <code>UNKNOWN</code> indicates that the alignment format is
     * unknown.
     */
    public static final int UNKNOWN = 100;

    /**
     * <code>RAW</code> indicates that the alignment format is raw
     * (symbols only).
     */
    public static final int RAW = 101;

    /**
     * <code>FASTA</code> indicates that the alignment format is
     * Fasta.
     */
    public static final int FASTA = 102;

    /**
     * <code>CLUSTAL</code> indicates that the alignment format is
     * Clustal.
     */
    public static final int CLUSTAL = 103;

    /**
     * <code>MSF</code> indicates that the alignment format is MSF.
     */
    public static final int MSF = 104;

    /**
     * <code>RAW_DNA</code> premade RAW | DNA.
     */
    public static final int RAW_DNA = RAW | SeqIOConstants.DNA;

    /**
     * <code>RAW_RNA</code> premade RAW | RNA.
     */
    public static final int RAW_RNA = RAW | SeqIOConstants.RNA;

    /**
     * <code>RAW_AA</code> premade RAW | AA.
     */
    public static final int RAW_AA = RAW | SeqIOConstants.AA;

    /**
     * <code>FASTA_DNA</code> premade FASTA | DNA;
     */
    public static final int FASTA_DNA = FASTA | SeqIOConstants.DNA;

    /**
     * <code>FASTA_RNA</code> premade FASTA | RNA;
     */
    public static final int FASTA_RNA = FASTA | SeqIOConstants.RNA;

    /**
     * <code>FASTA_AA</code> premade FASTA | AA;
     */
    public static final int FASTA_AA = FASTA | SeqIOConstants.AA;

    /**
     * <code>CLUSTAL_DNA</code> premade CLUSTAL | DNA;
     */
    public static final int CLUSTAL_DNA = CLUSTAL | SeqIOConstants.DNA;

    /**
     * <code>CLUSTAL_RNA</code> premade CLUSTAL | RNA;
     */
    public static final int CLUSTAL_RNA = CLUSTAL | SeqIOConstants.RNA;

    /**
     * <code>CLUSTAL_AA</code> premade CLUSTAL | AA;
     */
    public static final int CLUSTAL_AA  = CLUSTAL | SeqIOConstants.AA;

    /**
     * <code>MSF_DNA</code> premade MSF | DNA;
     */
    public static final int MSF_DNA = MSF | SeqIOConstants.DNA;

    /**
     * <code>MSF_DNA</code> premade MSF | RNA;
     */
    public static final int MSF_RNA = MSF | SeqIOConstants.RNA;

    /**
     * <code>MSF_AA</code> premade MSF | AA;
     */
    public static final int MSF_AA = MSF | SeqIOConstants.AA;
}
