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

package org.biojava.bio.program.gff;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.utils.ParserException;

/**
 * Interface which captures any errors which occur when parsing
 * a GFF stream.  Providing a custom implementation of this
 * interface allows intelligent recovery from errors when
 * parsing GFF.
 *
 * <p>
 * Each of these methods has three options:
 *
 * <ul>
 * <li>Throw a ParserException.  This need only contain a 
 * detail message, the parser will fill in other fields.
 * parsing will be aborted.</li>
 * <li>Throw an IgnoreRecordException.  This line of the GFF
 * file will be ignored, but parsing will not be aborted</li>
 * <li>Return a value for the field.
 * </ul>
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public interface GFFErrorHandler {
    /** 
     * The `start' field of the GFF entry was not a valid value.
     *
     * @param token The start token found.
     * @return A parsed value, if this is possible
     * @throws ParserException If parsing should be aborted
     * @throws IgnoreRecordException If this record should be silently skipped.
     */

    public int invalidStart(String token) 
        throws ParserException, IgnoreRecordException;

    /** 
     * The `end' field of the GFF entry was not a valid value.
     *
     * @param token The end token found.
     * @return A parsed value, if this is possible
     * @throws ParserException If parsing should be aborted
     * @throws IgnoreRecordException If this record should be silently skipped.
     */

    public int invalidEnd(String token) 
        throws ParserException, IgnoreRecordException;

    /** 
     * The `score' field of the GFF entry was not a valid value.
     *
     * @param token The score token found.
     * @return A parsed value, if this is possible
     * @throws ParserException If parsing should be aborted
     * @throws IgnoreRecordException If this record should be silently skipped.
     */

    public double invalidScore(String token) 
        throws ParserException, IgnoreRecordException;

    /** 
     * The `frame' field of the GFF entry was not a valid value.
     *
     * @param token The frame token found.
     * @return A parsed value, if this is possible
     * @throws ParserException If parsing should be aborted
     * @throws IgnoreRecordException If this record should be silently skipped.
     */

    public int invalidFrame(String token) 
        throws ParserException, IgnoreRecordException;

    /** 
     * The `strand' field of the GFF entry was not a valid value.
     *
     * @param token The strand token found.
     * @return A parsed value, if this is possible
     * @throws ParserException If parsing should be aborted
     * @throws IgnoreRecordException If this record should be silently skipped.
     */

    public StrandedFeature.Strand invalidStrand(String token) 
        throws ParserException, IgnoreRecordException;

    public static final GFFErrorHandler ABORT_PARSING = new AbortErrorHandler();

    static class AbortErrorHandler implements GFFErrorHandler {
	public int invalidStart(String token) throws ParserException {
	    throw new ParserException("Invalid start token");
	}

	public int invalidEnd(String token) throws ParserException {
	    throw new ParserException("Invalid end token");
	}

	public double invalidScore(String token) throws ParserException {
	    throw new ParserException("Invalid score token");
	}

	public int invalidFrame(String token) throws ParserException {
	    throw new ParserException("Invalid frame token");
	}

	public StrandedFeature.Strand invalidStrand(String token) throws ParserException {
	    throw new ParserException("Invalid strand token");
	}
    }

    public static final GFFErrorHandler SKIP_RECORD = new SkipRecordErrorHandler();
    
    static class SkipRecordErrorHandler implements GFFErrorHandler {
	public int invalidStart(String token) throws IgnoreRecordException {
	    throw new IgnoreRecordException();
	}

	public int invalidEnd(String token) throws IgnoreRecordException {
	    throw new IgnoreRecordException();
	}
	
	public double invalidScore(String token) throws IgnoreRecordException {
	    throw new IgnoreRecordException();
	}

	public int invalidFrame(String token) throws IgnoreRecordException {
	    throw new IgnoreRecordException();
	}

	public StrandedFeature.Strand invalidStrand(String token) throws IgnoreRecordException {
	    throw new IgnoreRecordException();
	}
    }
}



