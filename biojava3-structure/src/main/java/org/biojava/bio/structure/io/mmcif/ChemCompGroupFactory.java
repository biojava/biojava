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
 * Created on May 23, 2010
 *
 */
package org.biojava.bio.structure.io.mmcif;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.AminoAcidImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.NucleotideImpl;

import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava3.core.util.SoftHashMap;


public class ChemCompGroupFactory {

	static ChemCompProvider chemCompProvider = new DownloadChemCompProvider();

	static SoftHashMap<String, ChemComp> cache = new SoftHashMap<String, ChemComp>(0);

	public static ChemComp getChemComp(String recordName){

		recordName = recordName.toUpperCase().trim();

		// we are using the cache, to avoid hitting the file system too often.
		ChemComp cc = cache.get(recordName);
		if ( cc != null)
			return cc;

		// not cached, get the chem comp from the provider 
		cc = chemCompProvider.getChemComp(recordName);
		
		cache.put(recordName, cc);
		return cc;
	}

	public static void setChemCompProvider(ChemCompProvider provider){
		chemCompProvider = provider;
	}

	public static ChemCompProvider getChemCompProvider(){
		return chemCompProvider;
	}

	public static Group getGroupFromChemCompDictionary(String recordName) {

		// make sure we work with upper case records		
		recordName = recordName.toUpperCase().trim();

		Group g = null;

		try {

			ChemComp cc =  getChemComp(recordName);

			if ( cc == null)
				return null;

			if ( PolymerType.PROTEIN_ONLY.contains( cc.getPolymerType() ) ){
				AminoAcid aa = new AminoAcidImpl();

				String one_letter = cc.getOne_letter_code();
				if ( one_letter == null || one_letter.equals("X") || one_letter.equals("?") || one_letter.length()==0){
					String parent = cc.getMon_nstd_parent_comp_id();
					if ( parent != null && parent.length() == 3){
						String parentid = cc.getMon_nstd_parent_comp_id() ;
						ChemComp parentCC = getChemComp(parentid);
						one_letter = parentCC.getOne_letter_code();
					}
				}
				
				if ( one_letter == null || one_letter.length()==0 || one_letter.equals("?")) {
					// e.g. problem with PRR, which probably should have a parent of ALA, but as of 20110127 does not.
					System.err.println(" Problem with chemical component: " + recordName + "  Did not find one letter code!");
					aa.setAminoType('X');

				} else  {
					aa.setAminoType(one_letter.charAt(0));
				}


				g = aa;
			} else if ( PolymerType.POLYNUCLEOTIDE_ONLY.contains(cc.getPolymerType())) {
				NucleotideImpl nuc = new NucleotideImpl();

				g = nuc;


			} else {

				g = new HetatomImpl();
			}

			g.setChemComp(cc);

		} catch (Exception e){
			e.printStackTrace();
		}
		return g;
	}


	public  static String getOneLetterCode(ChemComp cc){
		String oneLetter = cc.getOne_letter_code();
		if ( oneLetter == null || oneLetter.equals("X") || oneLetter.equals("?")) {
			String parentId = cc.getMon_nstd_parent_comp_id() ;
			if ( parentId == null)
				return oneLetter;
			ChemComp parentCC = ChemCompGroupFactory.getChemComp(parentId);
			if ( parentCC == null)
				return oneLetter;
			oneLetter = parentCC.getOne_letter_code();
		}
		return oneLetter;
	}




}
