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
 * Created on 2012-02-11
 *
 */

package org.biojava3.ws.alignment.qblast;

/**
 * Alignment request parameters accepted by QBlast service.<br/>
 * Not all are mandatory. Certain parameters only work with a subset of other parameters in the list.
 * <p/>
 * Taken from <a href=http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node9.html>Blast URL API</a>
 * 
 * @author Gediminas Rimsa
 */
public enum BlastAlignmentParameterEnum {

	CMD,
	TOOL,
	EMAIL,

	ALIGNMENTS,
	BLAST_PROGRAM,
	CDD_SEARCH,
	COMPOSITION_BASED_STATISTICS,
	DATABASE_PREFIX,
	DATABASE,
	DB_GENETIC_CODE,
	DESCRIPTIONS,
	ENTREZ_QUERY,
	EXPECT,
	FILTER,
	FIRST_QUERY_NUM,
	GAPCOSTS,
	GENETIC_CODE,
	HITLIST_SIZE,
	I_THRESH,
	LCASE_MASK,
	MATCH_SCORES,
	MATRIX_NAME,
	MAX_NUM_SEQ,
	MEGABLAST,
	NUM_OVERVIEW,
	OTHER_ADVANCED,
	PERC_IDENT,
	PHI_PATTERN,
	PROGRAM,
	PSSM,
	QUERY,
	QUERY_BELIEVE_DEFLINE,
	QUERY_FROM,
	QUERY_TO,
	REPEATS,
	SHORT_QUERY_ADJUST,
	SEARCHSP_EFF,
	TEMPLATE_LENTH,
	TEMPLATE_TYPE,
	THRESHOLD,
	TWO_HITS,
	WORD_SIZE,
	WWW_BLAST_TYPE

}
