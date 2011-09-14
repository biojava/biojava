package org.biojava.bio.structure;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Scanner;

import junit.framework.TestCase;

public abstract class StringManipulationTestsHelper extends TestCase {

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
				sb.append(line + "\n");
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
