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
package org.biojava.nbio.structure.io;

//import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;

import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.junit.Test;

public class TestMmCIFSpecialCases {

	/**
	 * This tests for cases where dots appear in integer fields.
	 * Unusual but it happens in some PDB entries like 1s32
	 * See issue https://github.com/biojava/biojava/issues/368
	 * @throws IOException
	 */
	@Test
	public void testDotsInIntFields() throws IOException {

		// taken from 1s32
		String mmcifStr =
				"data_\n" +
				"loop_\n" +
				"_struct_ref_seq_dif.align_id\n" +
				"_struct_ref_seq_dif.pdbx_pdb_id_code\n"+
				"_struct_ref_seq_dif.mon_id\n"+
				"_struct_ref_seq_dif.pdbx_pdb_strand_id\n"+
				"_struct_ref_seq_dif.seq_num\n"+ // integer field that contains '.'
				"_struct_ref_seq_dif.pdbx_seq_db_name\n"+
				"_struct_ref_seq_dif.pdbx_seq_db_accession_code\n"+
				"_struct_ref_seq_dif.db_mon_id\n"+
				"_struct_ref_seq_dif.pdbx_seq_db_seq_num\n"+
				"_struct_ref_seq_dif.details\n"+
				"_struct_ref_seq_dif.pdbx_auth_seq_num\n"+
				"_struct_ref_seq_dif.pdbx_pdb_ins_code\n"+
				"_struct_ref_seq_dif.pdbx_ordinal\n"+
				"1 1S32 . A . GB  30268544 MET 1 'INTIATING METHIONINE' ? ? 1\n"+
				"2 1S32 . E . GB  30268544 MET 1 'INTIATING METHIONINE' ? ? 2\n"+
				"3 1S32 . B . UNP P02304   MET 0 'INTIATING METHIONINE' ? ? 3\n"+
				"4 1S32 . F . UNP P02304   MET 0 'INTIATING METHIONINE' ? ? 4\n"+
				"5 1S32 . C . GB  30268540 MET 1 'INTIATING METHIONINE' ? ? 5\n"+
				"6 1S32 . G . GB  30268540 MET 1 'INTIATING METHIONINE' ? ? 6\n"+
				"7 1S32 . D . GB  30268542 MET 1 'INTIATING METHIONINE' ? ? 7\n"+
				"8 1S32 . H . GB  30268542 MET 1 'INTIATING METHIONINE' ? ? 8" ;

		SimpleMMcifParser parser = new SimpleMMcifParser();

		BufferedReader buf = new BufferedReader(new StringReader(mmcifStr));

		parser.parse(buf);

		buf.close();

		// nothing to assert, the test just makes sure it doesn't throw an exception


	}

}
