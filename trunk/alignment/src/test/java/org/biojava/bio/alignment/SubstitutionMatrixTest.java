/*
 *                  BioJava development code
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
 * Created on Aug 23, 2007
 *
 */

package org.biojava.bio.alignment;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import junit.framework.TestCase;

import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;

public class SubstitutionMatrixTest extends TestCase {



	public void testParseSubstitutionMatrix(){

		InputStream inStream = this.getClass().getResourceAsStream("/blosum62.mat");
        assertNotNull(inStream);


        try {
            FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
            SymbolTokenization symtok = alphabet.getTokenization("token");
        	//String file = readMatrix(inStream);
        	//SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet,file,"blosum 62");
        	SubstitutionMatrix matrix = SubstitutionMatrix.getSubstitutionMatrix(
        			new BufferedReader(new InputStreamReader(inStream)));
        	//matrix.printMatrix();

        	Symbol A = symtok.parseToken("A");
        	Symbol W = symtok.parseToken("W");
        	Symbol D = symtok.parseToken("D");


        	assertEquals(matrix.getValueAt(A, A), 4);
        	assertEquals(matrix.getValueAt(W, D),-4);
        } catch (Exception e){
        	fail(e.getMessage());
        }


	}


	private String readMatrix(InputStream stream) throws IOException{
		String newline = System.getProperty("line.separator");
		BufferedReader reader = new BufferedReader(new InputStreamReader( stream));
		StringBuffer file = new StringBuffer();
		while (reader.ready()){
			file.append(reader.readLine() );
			file.append(newline);
		}

		return file.toString();
	}
}
