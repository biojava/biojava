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
 * Created on 2011-11-20
 *
 */

package org.biojava3.ws.alignment.qblast2.enums;

/**
 * Some constants representing available blast parameters (the list is not
 * complete!)
 * <p>
 * To easily use them in your application add a static import:
 * <p>
 * {@code import static full.package.name.BlastParameter.*;}
 * 
 * @author Gediminas Rimsa
 */
public abstract class BlastParameter {
	// Generic parameters
	public static final String CMD = "CMD";
	public static final String RID = "RID";
	public static final String FORMAT_TYPE = "FORMAT_TYPE";

	// Alignment request parameters
	public static final String QUERY = "QUERY";
	public static final String PROGRAM = "PROGRAM";
	public static final String MEGABLAST = "MEGABLAST";
	public static final String DATABASE = "DATABASE";
	public static final String EXPECT = "EXPECT";
	public static final String WORD_SIZE = "WORD_SIZE";
	public static final String GAPCOSTS = "GAPCOSTS";
	public static final String MATRIX_NAME = "MATRIX_NAME";

	
	public static final String QUERY_FROM = "QUERY_FROM";
	public static final String QUERY_TO = "QUERY_TO";
	public static final String OTHER_ADVANCED = "OTHER_ADVANCED";
	public static final String ENTREZ_QUERY = "ENTREZ_QUERY";
	public static final String TOOL = "TOOL";
	public static final String EMAIL = "EMAIL";

	// Output parameters
	public static final String ALIGNMENT_VIEW = "ALIGNMENT_VIEW";
	public static final String DESCRIPTIONS = "DESCRIPTIONS";
	public static final String ALIGNMENTS = "ALIGNMENTS";

}
