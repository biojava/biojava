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
 * created at Mar 4, 2008
 */
package org.biojava.nbio.structure.io.mmcif;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.MMCIFFileReader;
import org.biojava.nbio.structure.io.StructureIOFile;
import org.biojava.nbio.structure.io.mmcif.model.AtomSite;
import org.biojava.nbio.structure.io.mmcif.model.AtomSites;
import org.biojava.nbio.structure.io.mmcif.model.AuditAuthor;
import org.biojava.nbio.structure.io.mmcif.model.CIFLabel;
import org.biojava.nbio.structure.io.mmcif.model.Cell;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompAtom;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompBond;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompDescriptor;
import org.biojava.nbio.structure.io.mmcif.model.DatabasePDBremark;
import org.biojava.nbio.structure.io.mmcif.model.DatabasePDBrev;
import org.biojava.nbio.structure.io.mmcif.model.DatabasePdbrevRecord;
import org.biojava.nbio.structure.io.mmcif.model.Entity;
import org.biojava.nbio.structure.io.mmcif.model.EntityPoly;
import org.biojava.nbio.structure.io.mmcif.model.EntityPolySeq;
import org.biojava.nbio.structure.io.mmcif.model.EntitySrcGen;
import org.biojava.nbio.structure.io.mmcif.model.EntitySrcNat;
import org.biojava.nbio.structure.io.mmcif.model.EntitySrcSyn;
import org.biojava.nbio.structure.io.mmcif.model.Exptl;
import org.biojava.nbio.structure.io.mmcif.model.IgnoreField;
import org.biojava.nbio.structure.io.mmcif.model.PdbxAuditRevisionHistory;
import org.biojava.nbio.structure.io.mmcif.model.PdbxChemCompDescriptor;
import org.biojava.nbio.structure.io.mmcif.model.PdbxChemCompIdentifier;
import org.biojava.nbio.structure.io.mmcif.model.PdbxDatabaseStatus;
import org.biojava.nbio.structure.io.mmcif.model.PdbxEntityNonPoly;
import org.biojava.nbio.structure.io.mmcif.model.PdbxNonPolyScheme;
import org.biojava.nbio.structure.io.mmcif.model.PdbxPolySeqScheme;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructOperList;
import org.biojava.nbio.structure.io.mmcif.model.Refine;
import org.biojava.nbio.structure.io.mmcif.model.Struct;
import org.biojava.nbio.structure.io.mmcif.model.StructAsym;
import org.biojava.nbio.structure.io.mmcif.model.StructConn;
import org.biojava.nbio.structure.io.mmcif.model.StructKeywords;
import org.biojava.nbio.structure.io.mmcif.model.StructNcsOper;
import org.biojava.nbio.structure.io.mmcif.model.StructRef;
import org.biojava.nbio.structure.io.mmcif.model.StructRefSeq;
import org.biojava.nbio.structure.io.mmcif.model.StructRefSeqDif;
import org.biojava.nbio.structure.io.mmcif.model.StructSite;
import org.biojava.nbio.structure.io.mmcif.model.StructSiteGen;
import org.biojava.nbio.structure.io.mmcif.model.Symmetry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A simple mmCif file parser
 *
 *
 * Usage:
 * <pre>
String file = "path/to/mmcif/file";
StructureIOFile pdbreader = new MMCIFFileReader();

Structure s = pdbreader.getStructure(file);
System.out.println(s);

// you can convert it to a PDB file...
System.out.println(s.toPDB());

 * </pre>
 * For more documentation see <a href="http://biojava.org/wiki/BioJava:CookBook#Protein_Structure">http://biojava.org/wiki/BioJava:CookBook#Protein_Structure</a>.
 *
 * @author Andreas Prlic
 * @author Jose Duarte
 * @since 1.7
 */
public class SimpleMMcifParser implements MMcifParser {



	/**
	 * The header appearing at the beginning of a mmCIF file.
	 * A "block code" can be added to it of no more than 32 chars.
	 * See http://www.iucr.org/__data/assets/pdf_file/0019/22618/cifguide.pdf
	 */
	public static final String MMCIF_TOP_HEADER = "data_";

	public static final String COMMENT_CHAR = "#";
	public static final String LOOP_START = "loop_";
	public static final String FIELD_LINE = "_";

	// the following are the 3 valid quoting characters in CIF
	/**
	 * Quoting character '
	 */
	private static final char S1 = '\'';

	/**
	 * Quoting character "
	 */
	private static final char S2 = '\"';

	/**
	 * Quoting character ; (multi-line quoting)
	 */
	public static final String STRING_LIMIT = ";";


	private List<MMcifConsumer> consumers ;

	private Struct struct ;

	private static final Logger logger = LoggerFactory.getLogger(SimpleMMcifParser.class);

	public SimpleMMcifParser(){
		consumers = new ArrayList<MMcifConsumer>();
		struct = null;
	}

	@Override
	public void addMMcifConsumer(MMcifConsumer consumer) {
		consumers.add(consumer);

	}

	@Override
	public void clearConsumers() {
		consumers.clear();

	}

	@Override
	public void removeMMcifConsumer(MMcifConsumer consumer) {
		consumers.remove(consumer);
	}

	public static void main(String[] args){
		String file = "/Users/andreas/WORK/PDB/mmCif/a9/1a9n.cif.gz";
		//String file = "/Users/andreas/WORK/PDB/MMCIF/1gav.mmcif";
		//String file = "/Users/andreas/WORK/PDB/MMCIF/100d.cif";
		//String file = "/Users/andreas/WORK/PDB/MMCIF/1a4a.mmcif";
		System.out.println("parsing " + file);

		StructureIOFile pdbreader = new MMCIFFileReader();
		try {
			Structure s = pdbreader.getStructure(file);
			System.out.println(s);
			// convert it to a PDB file...
			System.out.println(s.toPDB());
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	@Override
	public void parse(InputStream inStream) throws IOException {
		parse(new BufferedReader(new InputStreamReader(inStream)));

	}

	@Override
	public void parse(BufferedReader buf)
			throws IOException {

		triggerDocumentStart();


		// init container objects...
		struct = new Struct();
		String line = null;

		boolean inLoop = false;
		boolean inLoopData = false;


		List<String> loopFields = new ArrayList<String>();
		List<String> lineData   = new ArrayList<String>();
		Set<String> loopWarnings = new HashSet<String>(); // used only to reduce logging statements

		String category = null;
		
		boolean foundHeader = false;

		while ( (line = buf.readLine ()) != null ){

			if (line.isEmpty() || line.startsWith(COMMENT_CHAR)) continue;

			if (!foundHeader) {
				// the first non-comment line is a data_PDBCODE line, test if this looks like a mmcif file
				if (line.startsWith(MMCIF_TOP_HEADER)){
					foundHeader = true;
					continue;
				} else {
					triggerDocumentEnd();
					throw new IOException("This does not look like a valid mmCIF file! The first line should start with 'data_', but is: '" + line+"'");
				}
			}

			logger.debug(inLoop + " " + line);

			if (line.startsWith(MMCIF_TOP_HEADER)){
				// either first line in file, or beginning of new section (data block in CIF parlance)
				if ( inLoop) {
					//System.out.println("new data and in loop: " + line);
					inLoop = false;
					inLoopData = false;
					lineData.clear();
					loopFields.clear();
				}

			}


			if ( inLoop) {


				if ( line.startsWith(LOOP_START)){
					loopFields.clear();
					inLoop = true;
					inLoopData = false;
					continue;
				}

				if ( line.matches("\\s*"+FIELD_LINE+"\\w+.*")) {

					if (inLoopData && line.startsWith(FIELD_LINE)) {
						logger.debug("Found a field line after reading loop data. Toggling to inLoop=false");
						inLoop = false;
						inLoopData = false;
						loopFields.clear();


						// a boring normal line
						List<String> data = processLine(line, buf, 2);

						if ( data.size() < 1){
							// this can happen if empty lines at end of file
							lineData.clear();
							continue;
						}
						String key = data.get(0);
						int pos = key.indexOf(".");
						if ( pos < 0 ) {
							// looks like a chem_comp file
							// line should start with data, otherwise something is wrong!
							if (! line.startsWith(MMCIF_TOP_HEADER)){
								logger.warn("This does not look like a valid mmCIF file! The first line should start with 'data_', but is '" + line+"'");
								triggerDocumentEnd();
								return;
							}
							// ignore the first line...
							category=null;
							lineData.clear();
							continue;
						}
						category = key.substring(0,pos);
						String value = data.get(1);
						loopFields.add(key.substring(pos+1,key.length()));
						lineData.add(value);

						logger.debug("Found data for category {}: {}", key, value);
						continue;
					}

					// found another field.
					String txt = line.trim();
					if ( txt.indexOf('.') > -1){

						String[] spl = txt.split("\\.");
						category = spl[0];
						String attribute = spl[1];
						loopFields.add(attribute);
						logger.debug("Found category: {}, attribute: {}",category, attribute);
						if ( spl.length > 2){
							logger.warn("Found nested attribute in {}, not supported yet!",txt);
						}

					} else {
						category = txt;
						logger.debug("Found category without attribute: {}",category);
					}


				} else {

					// in loop and we found a data line
					lineData = processLine(line, buf, loopFields.size());
					logger.debug("Found a loop data line with {} data fields", lineData.size());
					logger.debug("Data fields: {}", lineData.toString());
					if ( lineData.size() != loopFields.size()){
						logger.warn("Expected {} data fields, but found {} in line: {}",loopFields.size(),lineData.size(),line);

					}

					endLineChecks(category, loopFields, lineData, loopWarnings);

					lineData.clear();

					inLoopData = true;
				}

			} else {
				// not in loop

				if ( line.startsWith(LOOP_START)){
					if ( category != null)
						endLineChecks(category, loopFields, lineData, loopWarnings);

					resetBuffers(loopFields, lineData, loopWarnings);
					category = null;
					inLoop = true;
					inLoopData = false;
					logger.debug("Detected LOOP_START: '{}'. Toggling to inLoop=true", LOOP_START);
					continue;
				} else {
					logger.debug("Normal line ");
					inLoop = false;

					// a boring normal line
					List<String> data = processLine(line, buf, 2);

					if ( data.size() < 1){
						// this can happen if empty lines at end of file
						lineData.clear();
						continue;
					}
					String key = data.get(0);
					int pos = key.indexOf(".");
					if ( pos < 0 ) {
						// looks like a chem_comp file
						// line should start with data, otherwise something is wrong!
						if (! line.startsWith(MMCIF_TOP_HEADER)){
							logger.warn("This does not look like a valid mmCIF file! The first line should start with 'data_', but is '" + line+"'");
							triggerDocumentEnd();
							return;
						}
						// ignore the first line...
						category=null;
						lineData.clear();
						continue;
					}

					if (category!=null && !key.substring(0,pos).equals(category)) {
						// we've changed category: need to flush the previous one
						endLineChecks(category, loopFields, lineData, loopWarnings);
						resetBuffers(loopFields, lineData, loopWarnings);
					}

					category = key.substring(0,pos);

					String value = data.get(1);
					loopFields.add(key.substring(pos+1,key.length()));
					lineData.add(value);

					logger.debug("Found data for category {}: {}", key, value);

				}
			}
		}

		if (category!=null && lineData.size()>0 && lineData.size()==loopFields.size()) {
			// the last category in the file will still be missing, we add it now
			endLineChecks(category, loopFields, lineData, loopWarnings);
			resetBuffers(loopFields, lineData, loopWarnings);
		}

		if (struct != null){
			triggerStructData(struct);
		}

		triggerDocumentEnd();

	}

	private void resetBuffers(List<String> loopFields, List<String> lineData, Set<String> loopWarnings) {
		loopFields.clear();
		lineData.clear();
		loopWarnings.clear();
	}

	private List<String> processSingleLine(String line){

		List<String> data = new ArrayList<String>();

		if ( line.trim().length() == 0){
			return data;
		}

		if ( line.trim().length() == 1){
			if ( line.startsWith(STRING_LIMIT))
				return data;
		}
		boolean inString = false; // semicolon (;) quoting
		boolean inS1     = false; // single quote (') quoting
		boolean inS2     = false; // double quote (") quoting
		String word 	 = "";

		for (int i=0; i< line.length(); i++ ){

			Character c = line.charAt(i);

			Character nextC = null;
			if (i < line.length() - 1)
				nextC = line.charAt(i+1);

			Character prevC = null;
			if (i>0)
				prevC = line.charAt(i-1);

			if  (c == ' ') {

				if ( ! inString){
					if ( ! word.equals(""))
						data.add(word.trim());
					word = "";
				} else {
					// we are in a string, add the space
					word += c;
				}

			} else if (c == S1 )  {

				if ( inString){

					boolean wordEnd = false;
					if (! inS2) {
						if (nextC==null || Character.isWhitespace(nextC)){
							i++;
							wordEnd = true;
						}
					}


					if ( wordEnd ) {

						// at end of string
						if ( ! word.equals(""))
							data.add(word.trim());
						word     = "";
						inString = false;
						inS1     = false;
					} else {
						word += c;
					}

				} else if (prevC==null || prevC==' ') {
					// the beginning of a new string
					inString = true;
					inS1     = true;
				} else {
					word += c;
				}
			} else if ( c == S2 ){
				if ( inString){

					boolean wordEnd = false;
					if (! inS1) {
						if (nextC==null || Character.isWhitespace(nextC)){
							i++;
							wordEnd = true;
						}
					}

					if ( wordEnd ) {

						// at end of string
						if ( ! word.equals(""))
							data.add(word.trim());
						word     = "";
						inString = false;
						inS2     = false;
					} else {
						word += c;
					}
				}  else if (prevC==null || prevC==' ') {
					// the beginning of a new string
					inString = true;
					inS2     = true;
				} else {
					word += c;
				}
			} else {
				word += c;
			}

		}
		if ( ! word.trim().equals(""))
			data.add(word);


		return data;

	}

	/**
	 * Get the content of a cif entry
	 *
	 * @param line
	 * @param buf
	 * @return
	 */
	private List<String> processLine(String line,
			BufferedReader buf,
			int fieldLength)
					throws IOException{

		//System.out.println("XX processLine " + fieldLength + " " + line);
		// go through the line and process each character
		List<String> lineData = new ArrayList<String>();

		boolean inString = false;

		StringBuilder bigWord = null;

		while ( true ){

			if ( line.startsWith(STRING_LIMIT)){
				if (! inString){

					inString = true;
					if ( line.length() > 1)
						bigWord = new StringBuilder(line.substring(1));
					else
						bigWord = new StringBuilder("");


				} else {
					// the end of a word
					lineData.add(bigWord.toString());
					bigWord = null;
					inString = false;

				}
			} else {
				if ( inString )
					bigWord.append(line);
				else {

					List<String> dat = processSingleLine(line);

					for (String d : dat){
						lineData.add(d);
					}
				}
			}

			//System.out.println("in process line : " + lineData.size() + " " + fieldLength);

			if ( lineData.size() > fieldLength){

				logger.warn("wrong data length ("+lineData.size()+
						") should be ("+fieldLength+") at line " + line + " got lineData: " + lineData);
				return lineData;
			}

			if ( lineData.size() == fieldLength)
				return lineData;


			line = buf.readLine();
			if ( line == null)
				break;
		}
		return lineData;

	}



	private void endLineChecks(String category,List<String> loopFields, List<String> lineData, Set<String> loopWarnings ) throws IOException{

		logger.debug("Processing category {}, with fields: {}",category,loopFields.toString());
		//		System.out.println("parsed the following data: " +category + " fields: "+
		//				loopFields + " DATA: " +
		//				lineData);

		if ( loopFields.size() != lineData.size()){
			logger.warn("looks like we got a problem with nested string quote characters:");
			throw new IOException("data length ("+ lineData.size() +
					") != fields length ("+loopFields.size()+
					") category: " +category + " fields: "+
					loopFields + " DATA: " +
					lineData );
		}

		if ( category.equals("_entity")){

			Entity e =  (Entity) buildObject(
					Entity.class.getName(),
					loopFields,lineData, loopWarnings);
			triggerNewEntity(e);

		} else if (category.equals("_entity_poly")) {
			EntityPoly ep = (EntityPoly) buildObject(EntityPoly.class.getName(), loopFields, lineData, loopWarnings);
			triggerNewEntityPoly(ep);
			
		} else if ( category.equals("_struct")){

			struct =  (Struct) buildObject(
					Struct.class.getName(),
					loopFields, lineData, loopWarnings);

		} else if ( category.equals("_atom_site")){

			AtomSite a = (AtomSite) buildObject(
					AtomSite.class.getName(),
					loopFields, lineData, loopWarnings);
			triggerNewAtomSite(a);

		} else if ( category.equals("_database_PDB_rev")){
			DatabasePDBrev dbrev = (DatabasePDBrev) buildObject(
					DatabasePDBrev.class.getName(),
					loopFields, lineData, loopWarnings);

			triggerNewDatabasePDBrev(dbrev);

		} else if ( category.equals("_database_PDB_rev_record")) {
			DatabasePdbrevRecord dbrev = (DatabasePdbrevRecord) buildObject(
					DatabasePdbrevRecord.class.getName(),
					loopFields, lineData, loopWarnings);

			triggerNewDatabasePDBrevRecord(dbrev);
			
    // MMCIF version 5 dates  
		} else if ( category.equals("_pdbx_audit_revision_history")) {
			PdbxAuditRevisionHistory history = (PdbxAuditRevisionHistory) buildObject(
					PdbxAuditRevisionHistory.class.getName(),
					loopFields, lineData, loopWarnings);

			triggerNewPdbxAuditRevisionHistory(history);
    
    // MMCIF version 5 dates
		} else if ( category.equals("_pdbx_database_status")) {
			PdbxDatabaseStatus status = (PdbxDatabaseStatus) buildObject(
					PdbxDatabaseStatus.class.getName(),
					loopFields, lineData, loopWarnings);

			triggerNewPdbxDatabaseStatus(status);

		}else if (  category.equals("_database_PDB_remark")) {
			DatabasePDBremark remark = (DatabasePDBremark) buildObject(
					DatabasePDBremark.class.getName(),
					loopFields, lineData, loopWarnings);

			triggerNewDatabasePDBremark(remark);

		} else if ( category.equals("_exptl")){
			Exptl exptl  = (Exptl) buildObject(
					Exptl.class.getName(),
					loopFields,lineData, loopWarnings);

			triggerExptl(exptl);

		} else if ( category.equals("_cell")){
			Cell cell  = (Cell) buildObject(
					Cell.class.getName(),
					loopFields,lineData, loopWarnings);

			triggerNewCell(cell);

		} else if ( category.equals("_symmetry")){
			Symmetry symmetry  = (Symmetry) buildObject(
					Symmetry.class.getName(),
					loopFields,lineData, loopWarnings);

			triggerNewSymmetry(symmetry);
		} else if ( category.equals("_struct_ncs_oper")) {

			StructNcsOper sNcsOper = (StructNcsOper) buildObject(
					StructNcsOper.class.getName(), 
					loopFields, lineData, loopWarnings);
			triggerNewStructNcsOper(sNcsOper);
		} else if ( category.equals("_atom_sites")) {
			
			AtomSites atomSites = (AtomSites) buildObject(
					AtomSites.class.getName(),
					loopFields, lineData, loopWarnings);
			triggerNewAtomSites(atomSites);

		} else if ( category.equals("_struct_ref")){
			StructRef sref  = (StructRef) buildObject(
					StructRef.class.getName(),
					loopFields,lineData, loopWarnings);

			triggerNewStrucRef(sref);

		} else if ( category.equals("_struct_ref_seq")){
			StructRefSeq sref  = (StructRefSeq) buildObject(
					StructRefSeq.class.getName(),
					loopFields,lineData, loopWarnings);

			triggerNewStrucRefSeq(sref);
		} else if ( category.equals("_struct_ref_seq_dif")) {
			StructRefSeqDif sref = (StructRefSeqDif) buildObject(
					StructRefSeqDif.class.getName(),
					loopFields, lineData, loopWarnings);

			triggerNewStrucRefSeqDif(sref);
		} else if ( category.equals("_struct_site_gen")) {
			StructSiteGen sref = (StructSiteGen) buildObject(
					StructSiteGen.class.getName(),
					loopFields, lineData, loopWarnings);

			triggerNewStructSiteGen(sref);
		} else if ( category.equals("_struct_site")) {
			StructSite sref = (StructSite) buildObject(
					StructSite.class.getName(),
					loopFields, lineData, loopWarnings);
			triggerNewStructSite(sref);
		} else if ( category.equals("_entity_poly_seq")){
			EntityPolySeq exptl  = (EntityPolySeq) buildObject(
					EntityPolySeq.class.getName(),
					loopFields,lineData, loopWarnings);

			triggerNewEntityPolySeq(exptl);
		} else if ( category.equals("_entity_src_gen")){
			EntitySrcGen entitySrcGen = (EntitySrcGen) buildObject(
					EntitySrcGen.class.getName(),
					loopFields,lineData, loopWarnings);
			triggerNewEntitySrcGen(entitySrcGen);
		} else if ( category.equals("_entity_src_nat")){
			EntitySrcNat entitySrcNat = (EntitySrcNat) buildObject(
					EntitySrcNat.class.getName(),
					loopFields,lineData, loopWarnings);
			triggerNewEntitySrcNat(entitySrcNat);
		} else if ( category.equals("_pdbx_entity_src_syn")){
			EntitySrcSyn entitySrcSyn = (EntitySrcSyn) buildObject(
					EntitySrcSyn.class.getName(),
					loopFields,lineData, loopWarnings);
			triggerNewEntitySrcSyn(entitySrcSyn);
		} else if ( category.equals("_struct_asym")){
			StructAsym sasym  = (StructAsym) buildObject(
					StructAsym.class.getName(),
					loopFields,lineData, loopWarnings);

			triggerNewStructAsym(sasym);

		} else if ( category.equals("_pdbx_poly_seq_scheme")){
			PdbxPolySeqScheme ppss  = (PdbxPolySeqScheme) buildObject(
					PdbxPolySeqScheme.class.getName(),
					loopFields,lineData, loopWarnings);

			triggerNewPdbxPolySeqScheme(ppss);

		} else if ( category.equals("_pdbx_nonpoly_scheme")){
			PdbxNonPolyScheme ppss  = (PdbxNonPolyScheme) buildObject(
					PdbxNonPolyScheme.class.getName(),
					loopFields,lineData, loopWarnings);

			triggerNewPdbxNonPolyScheme(ppss);

		} else if ( category.equals("_pdbx_entity_nonpoly")){
			PdbxEntityNonPoly pen = (PdbxEntityNonPoly) buildObject(
					PdbxEntityNonPoly.class.getName(),
					loopFields,lineData, loopWarnings
					);
			triggerNewPdbxEntityNonPoly(pen);
		} else if ( category.equals("_struct_keywords")){
			StructKeywords kw = (StructKeywords)buildObject(
					StructKeywords.class.getName(),
					loopFields,lineData, loopWarnings
					);
			triggerNewStructKeywords(kw);
		} else if (category.equals("_refine")){
			Refine r = (Refine)buildObject(
					Refine.class.getName(),
					loopFields,lineData, loopWarnings
					);
			triggerNewRefine(r);
		} else if (category.equals("_chem_comp")){
			ChemComp c = (ChemComp)buildObject(
					ChemComp.class.getName(),
					loopFields, lineData, loopWarnings
					);
			triggerNewChemComp(c);
		} else if (category.equals("_audit_author")) {
			AuditAuthor aa = (AuditAuthor)buildObject(
					AuditAuthor.class.getName(),
					loopFields, lineData, loopWarnings);
			triggerNewAuditAuthor(aa);
		} else if (category.equals("_pdbx_chem_comp_descriptor")) {
			ChemCompDescriptor ccd = (ChemCompDescriptor) buildObject(
					ChemCompDescriptor.class.getName(),
					loopFields, lineData, loopWarnings);
			triggerNewChemCompDescriptor(ccd);
		} else if (category.equals("_pdbx_struct_oper_list")) {

			PdbxStructOperList structOper = (PdbxStructOperList) buildObject(
					PdbxStructOperList.class.getName(),
					loopFields, lineData, loopWarnings
					);
			triggerNewPdbxStructOper(structOper);

		} else if (category.equals("_pdbx_struct_assembly")) {
			PdbxStructAssembly sa = (PdbxStructAssembly) buildObject(
					PdbxStructAssembly.class.getName(),
					loopFields, lineData, loopWarnings);
			triggerNewPdbxStructAssembly(sa);

		} else if (category.equals("_pdbx_struct_assembly_gen")) {
			PdbxStructAssemblyGen sa = (PdbxStructAssemblyGen) buildObject(
					PdbxStructAssemblyGen.class.getName(),
					loopFields, lineData, loopWarnings);
			triggerNewPdbxStructAssemblyGen(sa);
		} else if ( category.equals("_chem_comp_atom")){
			ChemCompAtom atom = (ChemCompAtom)buildObject(
					ChemCompAtom.class.getName(),
					loopFields,lineData, loopWarnings);
			triggerNewChemCompAtom(atom);

		}else if ( category.equals("_chem_comp_bond")){
			ChemCompBond bond = (ChemCompBond)buildObject(
					ChemCompBond.class.getName(),
					loopFields,lineData, loopWarnings);
			triggerNewChemCompBond(bond);
		} else if ( category.equals("_pdbx_chem_comp_identifier")){
			PdbxChemCompIdentifier id = (PdbxChemCompIdentifier)buildObject(
					PdbxChemCompIdentifier.class.getName(),
					loopFields,lineData, loopWarnings);
			triggerNewPdbxChemCompIdentifier(id);
		} else if ( category.equals("_pdbx_chem_comp_descriptor")){
			PdbxChemCompDescriptor id = (PdbxChemCompDescriptor)buildObject(
					PdbxChemCompDescriptor.class.getName(),
					loopFields,lineData, loopWarnings);
			triggerNewPdbxChemCompDescriptor(id);
		} else if ( category.equals("_struct_conn")){
			StructConn id = (StructConn)buildObject(
					StructConn.class.getName(),
					loopFields,lineData, loopWarnings);
			triggerNewStructConn(id);

		} else {

			logger.debug("Using a generic bean for category {}",category);

			// trigger a generic bean that can deal with all missing data types...
			triggerGeneric(category,loopFields,lineData);
		}


	}


//	private PdbxStructOperList getPdbxStructOperList(List<String> loopFields,
//			List<String> lineData) {
//		PdbxStructOperList so = new PdbxStructOperList();
//
//		//System.out.println(loopFields);
//		//System.out.println(lineData);
//
//		String id = lineData.get(loopFields.indexOf("id"));
//		so.setId(id);
//		so.setType(lineData.get(loopFields.indexOf("type")));
//		Matrix matrix = new Matrix(3,3);
//		for (int i = 1 ; i <=3 ; i++){
//			for (int j =1 ; j <= 3 ; j++){
//				String max = String.format("matrix[%d][%d]",j,i);
//
//				String val = lineData.get(loopFields.indexOf(max));
//				Double d = Double.parseDouble(val);
//				matrix.set(j-1,i-1,d);
//				//				matrix.set(i-1,j-1,d);
//			}
//		}
//
//		double[] coords =new double[3];
//
//		for ( int i = 1; i <=3 ; i++){
//			String v = String.format("vector[%d]",i);
//			String val = lineData.get(loopFields.indexOf(v));
//			Double d = Double.parseDouble(val);
//			coords[i-1] = d;
//		}
//
//		so.setMatrix(matrix);
//		so.setVector(coords);
//
//
//
//		return so;
//	}

	public void triggerNewPdbxStructOper(PdbxStructOperList structOper) {
		for(MMcifConsumer c : consumers){
			c.newPdbxStructOperList(structOper);
		}

	}

	public void triggerNewStructNcsOper(StructNcsOper sNcsOper) {
		for(MMcifConsumer c : consumers){
			c.newStructNcsOper(sNcsOper);
		}

	}
	
	public void triggerNewAtomSites(AtomSites atomSites) {
		for(MMcifConsumer c : consumers){
			c.newAtomSites(atomSites);
		}
	}

	/**
	 * Populates a bean object from  the {@link org.biojava.nbio.structure.io.mmcif.model} package, 
	 * from the data read from a CIF file.
	 * It uses reflection to lookup the field and setter method names given the category 
	 * found in the CIF file. 
	 * <p>
	 * Due to limitations in variable names in java, not all fields can have names 
	 * exactly as defined in the CIF categories. In those cases the {@link CIFLabel} tag
	 * can be used in the field names to give the appropriate name that corresponds to the
	 * CIF category, which is the name that will be then looked up here.
	 * The {@link IgnoreField} tag can also be used to exclude fields from being looked up.
	 * @param className
	 * @param loopFields
	 * @param lineData
	 * @param warnings
	 * @return
	 */
	private Object buildObject(String className, List<String> loopFields, List<String> lineData, Set<String> warnings) {

		Object o = null;
		Class<?> c = null;

		try {
			// build up the Entity object from the line data...
			c = Class.forName(className);

			o = c.newInstance();

		} catch (InstantiationException|ClassNotFoundException|IllegalAccessException e){
			logger.error( "Error while constructing {}: {}", className, e.getMessage());
			return null;
		} 

		// these methods get the fields but also looking at the IgnoreField and CIFLabel annotations 
		Field[] fields = MMCIFFileTools.getFields(c);
		String[] names = MMCIFFileTools.getFieldNames(fields);

		// let's build a map of all methods so that we can look up the setter methods later
		Method[] methods = c.getMethods();

		Map<String,Method> methodMap = new HashMap<String, Method>();
		for (Method m : methods) {
			methodMap.put(m.getName(),m);
		}

		// and a map of all the fields so that we can lookup them up later
		Map<String, Field> names2fields = new HashMap<>();
		for (int i=0;i<fields.length;i++) {
			names2fields.put(names[i], fields[i]);
		}
		
		int pos = -1 ;
		for (String key: loopFields){
			pos++;

			String val = lineData.get(pos);
			
			// we first start looking up the field which can be annotated with a CIFLabel if they 
			// need alternative names (e.g. for field _symmetry.space_group_name_H-M, since hyphen is not allowed in var names in java)
			Field field = names2fields.get(key);
			
			if (field == null) {
				produceWarning(key, val, c, warnings);
				continue;
			}
			// now we need to find the corresponding setter
			// note that we can't use the field directly and then call Field.set() because many setters 
			// have more functionality than just setting the value (e.g. some setters in ChemComp)

			// building up the setter method name: need to upper case the first letter, leave the rest untouched
			String setterMethodName = "set" + field.getName().substring(0,1).toUpperCase() + field.getName().substring(1, field.getName().length());

			Method setter = methodMap.get(setterMethodName);
			
			if (setter==null) {
				produceWarning(key, val, c, warnings);
				continue;
			}
			
			

			// now we populate the object with the values by invoking the corresponding setter method,			
			// note that all of the mmCif container classes have only one argument (they are beans)
			Class<?>[] pType  = setter.getParameterTypes();
			

			try {
				if ( pType[0].getName().equals(Integer.class.getName())) {
					if ( val != null && ! val.equals("?") && !val.equals(".")) {

						Integer intVal = Integer.parseInt(val);
						setter.invoke(o, intVal);
						
					}
				} else {
					// default val is a String					
					setter.invoke(o, val);
				}
			} catch (IllegalAccessException|InvocationTargetException e) {
				logger.error("Could not invoke setter {} with value {} for class {}", setterMethodName, val, className);
			} 

		}

		return o;
	}
	
	private void produceWarning(String key, String val, Class<?> c, Set<String> warnings) {

		String warning = "Trying to set field " + key + " in "+ c.getName() +" found in file, but no corresponding field could be found in model class (value:" + val + ")";
		String warnkey = key+"-"+c.getName();
		// Suppress duplicate warnings or attempts to store empty data
		if( val.equals("?") || val.equals(".") || ( warnings != null && warnings.contains(warnkey)) ) {
			logger.debug(warning);
		} else {
			logger.info(warning);
		}

		if(warnings != null) {
			warnings.add(warnkey);
		}

	}

	public void triggerGeneric(String category, List<String> loopFields, List<String> lineData){
		for(MMcifConsumer c : consumers){
			c.newGenericData(category, loopFields, lineData);
		}
	}

	public void triggerNewEntity(Entity entity){
		for(MMcifConsumer c : consumers){
			c.newEntity(entity);
		}
	}
	
	public void triggerNewEntityPoly(EntityPoly entityPoly) {
		for(MMcifConsumer c : consumers){
			c.newEntityPoly(entityPoly);
		}		
	}

	public void triggerNewEntityPolySeq(EntityPolySeq epolseq){
		for(MMcifConsumer c : consumers){
			c.newEntityPolySeq(epolseq);
		}
	}
	public void triggerNewEntitySrcGen(EntitySrcGen entitySrcGen){
		for(MMcifConsumer c : consumers){
			c.newEntitySrcGen(entitySrcGen);
		}
	}
	public void triggerNewEntitySrcNat(EntitySrcNat entitySrcNat){
		for(MMcifConsumer c : consumers){
			c.newEntitySrcNat(entitySrcNat);
		}
	}
	public void triggerNewEntitySrcSyn(EntitySrcSyn entitySrcSyn){
		for(MMcifConsumer c : consumers){
			c.newEntitySrcSyn(entitySrcSyn);
		}
	}
	public void triggerNewChemComp(ChemComp cc){

		for(MMcifConsumer c : consumers){
			c.newChemComp(cc);
		}
	}
	public void triggerNewStructAsym(StructAsym sasym){
		for(MMcifConsumer c : consumers){
			c.newStructAsym(sasym);
		}
	}

	private void triggerStructData(Struct struct){
		for(MMcifConsumer c : consumers){
			c.setStruct(struct);
		}
	}

	private void triggerNewAtomSite(AtomSite atom){
		for(MMcifConsumer c : consumers){
			c.newAtomSite(atom);
		}
	}

	private void triggerNewAuditAuthor(AuditAuthor aa){
		for(MMcifConsumer c : consumers){
			c.newAuditAuthor(aa);
		}
	}
	
	private void triggerNewPdbxAuditRevisionHistory(PdbxAuditRevisionHistory history) {
		for(MMcifConsumer c : consumers){
			c.newPdbxAuditRevisionHistory(history);
		}
	}
	
	private void triggerNewPdbxDatabaseStatus(PdbxDatabaseStatus status) {
		for(MMcifConsumer c : consumers){
			c.newPdbxDatabaseStatus(status);
		}
	}
	
	private void triggerNewDatabasePDBrev(DatabasePDBrev dbrev){
		for(MMcifConsumer c : consumers){
			c.newDatabasePDBrev(dbrev);
		}
	}
	private void triggerNewDatabasePDBrevRecord(DatabasePdbrevRecord dbrev){
		for(MMcifConsumer c : consumers){
			c.newDatabasePDBrevRecord(dbrev);
		}
	}

	private void triggerNewDatabasePDBremark(DatabasePDBremark remark){
		for(MMcifConsumer c : consumers){
			c.newDatabasePDBremark(remark);
		}
	}

	private void triggerExptl(Exptl exptl){
		for(MMcifConsumer c : consumers){
			c.newExptl(exptl);
		}
	}

	private void triggerNewCell(Cell cell) {
		for(MMcifConsumer c : consumers){
			c.newCell(cell);
		}
	}

	private void triggerNewSymmetry(Symmetry symmetry) {
		for(MMcifConsumer c : consumers){
			c.newSymmetry(symmetry);
		}
	}

	private void triggerNewStrucRef(StructRef sref){
		for(MMcifConsumer c : consumers){
			c.newStructRef(sref);
		}
	}

	private void triggerNewStrucRefSeq(StructRefSeq sref){
		for(MMcifConsumer c : consumers){
			c.newStructRefSeq(sref);
		}
	}

	private void triggerNewStrucRefSeqDif(StructRefSeqDif sref){
		for(MMcifConsumer c : consumers){
			c.newStructRefSeqDif(sref);
		}
	}

	private void triggerNewPdbxPolySeqScheme(PdbxPolySeqScheme ppss){
		for(MMcifConsumer c : consumers){
			c.newPdbxPolySeqScheme(ppss);
		}
	}
	private void triggerNewPdbxNonPolyScheme(PdbxNonPolyScheme ppss){
		for(MMcifConsumer c : consumers){
			c.newPdbxNonPolyScheme(ppss);
		}
	}
	public void triggerNewPdbxEntityNonPoly(PdbxEntityNonPoly pen){
		for (MMcifConsumer c: consumers){
			c.newPdbxEntityNonPoly(pen);
		}
	}
	public void triggerNewStructKeywords(StructKeywords kw){
		for (MMcifConsumer c: consumers){
			c.newStructKeywords(kw);
		}
	}
	public void triggerNewRefine(Refine r){
		for (MMcifConsumer c: consumers){
			c.newRefine(r);
		}
	}
	public void triggerDocumentStart(){
		for(MMcifConsumer c : consumers){
			c.documentStart();
		}
	}
	public void triggerDocumentEnd(){
		for(MMcifConsumer c : consumers){
			c.documentEnd();
		}
	}
	public void triggerNewChemCompDescriptor(ChemCompDescriptor ccd) {
		for(MMcifConsumer c : consumers){
			c.newChemCompDescriptor(ccd);
		}
	}
	private void triggerNewPdbxStructAssembly(PdbxStructAssembly sa) {
		for(MMcifConsumer c : consumers){
			c.newPdbxStrucAssembly(sa);
		}
	}
	private void triggerNewPdbxStructAssemblyGen(PdbxStructAssemblyGen sa) {
		for(MMcifConsumer c : consumers){
			c.newPdbxStrucAssemblyGen(sa);
		}
	}

	private void triggerNewChemCompAtom(ChemCompAtom atom) {
		for(MMcifConsumer c : consumers){
			c.newChemCompAtom(atom);
		}
	}

	private void triggerNewChemCompBond(ChemCompBond bond) {
		for(MMcifConsumer c : consumers){
			c.newChemCompBond(bond);
		}
	}

	private void triggerNewPdbxChemCompIdentifier(PdbxChemCompIdentifier id) {
		for(MMcifConsumer c : consumers){
			c.newPdbxChemCompIndentifier(id);
		}
	}
	private void triggerNewPdbxChemCompDescriptor(PdbxChemCompDescriptor id) {
		for(MMcifConsumer c : consumers){
			c.newPdbxChemCompDescriptor(id);
		}
	}
	private void triggerNewStructConn(StructConn id) {
		for(MMcifConsumer c : consumers){
			c.newStructConn(id);
		}
	}
	private void triggerNewStructSiteGen(StructSiteGen id) {
		for (MMcifConsumer c : consumers) {
			c.newStructSiteGen(id);
		}
	}
	private void triggerNewStructSite(StructSite id) {
		for (MMcifConsumer c : consumers) {
			c.newStructSite(id);
		}
	}
}
