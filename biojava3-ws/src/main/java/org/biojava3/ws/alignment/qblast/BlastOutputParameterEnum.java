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
 * Output parameters accepted by QBlast service. <br/>
 * Not all are mandatory. Certain parameters only work with a subset of other parameters in the list.
 * <p/>
 * Taken from <a href=http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node9.html>Blast URL API</a>
 * 
 * @author Gediminas Rimsa
 */
public enum BlastOutputParameterEnum {

	ALIGNMENTS,
	ALIGNMENT_VIEW,
	DESCRIPTIONS,
	DATABASE_SORT,
	DISPLAY_SORT,
	EXPECT_HIGH,
	EXPECT_LOW,
	FORMAT_OBJECT,
	FORMAT_ENTREZ_QUERY,
	FORMAT_TYPE,
	HSP_SORT,
	I_THRESH,
	MASK_CHAR,
	MASK_COLOR,
	NCBI_GI,
	NOHEADER,
	NUM_OVERVIEW,
	RID,
	SHOW_CDS_FEATURE,
	SHOW_LINKOUT,
	SHOW_OVERVIEW

}
