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
 * Created on Jul 31, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.io.File;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;

import org.biojava3.protmod.parser.DefaultProteinModificationParser;

import junit.framework.TestCase;

public class ModifiedCompoundSerializationTest extends TestCase {

	public void testSerialization() throws StructureException, IOException, ClassNotFoundException {
		String pdbId = "3M6S";
		Structure struc = TmpAtomCache.cache.getStructure(pdbId);
		
		DefaultProteinModificationParser parser = new DefaultProteinModificationParser();
		parser.parse(struc, ProteinModification.allModifications());
		List<ModifiedCompound> mcs = parser.getIdentifiedModifiedCompound();
		
		String file = System.getProperty("java.io.tmpdir") + File.separatorChar + pdbId;
		FileOutputStream fos = new FileOutputStream(file);
		ObjectOutputStream oos = new ObjectOutputStream(fos);
		oos.writeObject(mcs);
		oos.close();
		
		FileInputStream fin = new FileInputStream(file);
		ObjectInputStream ois = new ObjectInputStream(fin); 
		mcs = (List<ModifiedCompound>) ois.readObject();
		ois.close();
		
		try {
			new File(file).delete();
		} catch (Exception e) {
			
		}
		
		//System.out.println(mcs);
	}
	
}
