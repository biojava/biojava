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
 * Created on Sep 14, 2011
 * Author: Amr AL-Hossary 
 *
 */
package org.biojava.bio.structure;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Scanner;

import junit.framework.TestCase;

/**A utility class for common string manipulation tasks.
 * <B>put here till it an be moved to some other more suitable package</B>
 * @author Amr AL-Hossary
 */
public abstract class StringManipulationTestsHelper extends TestCase {

	// we are using unix endline here, since this is used for testing XML and it is part of the XML recommendations:
	// http://www.w3.org/TR/REC-xml/#sec-line-ends
	private static final String UNIX_NEWLINE = "\n";
	
	private StringManipulationTestsHelper() {
//		to prevent instantiation
	}

	public static void assertEqualsIgnoreEndline(String str1, String str2) {
		Scanner scanner1 = new Scanner(str1);
		Scanner scanner2 = new Scanner(str2);
		String line1, line2;
		while (scanner1.hasNextLine()) {
			line1 = scanner1.nextLine();
			line2 = scanner2.nextLine();
			assertEquals(line1, line2);
		}
		if (scanner2.hasNextLine()) {
			//force fail
			assertEquals(str1, str2);
		}
	}

	public static void compareString(String t, String pdb){
		for (int i =0 ; i < t.length() ; i++){
			System.out.println("@"+i+"\t>"+t.charAt(i)+":"+ pdb.charAt(i)+"<\t"+Integer.toHexString(t.charAt(i))+":"+Integer.toHexString(pdb.charAt(i)));
			if ( Character.toUpperCase(t.charAt(i)) != Character.toUpperCase(pdb.charAt(i))){
				break;
			}
		}
	}

	public static String convertStreamToString(InputStream stream){
		BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
		StringBuilder sb = new StringBuilder();
	
		String line = null;
		try {
			while ((line = reader.readLine()) != null) {
			
				sb.append(line + UNIX_NEWLINE);
			}
		} catch (IOException e) {
			//e.printStackTrace();
		} finally {
			try {
				stream.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	
		return sb.toString();
	}
	
	

}
