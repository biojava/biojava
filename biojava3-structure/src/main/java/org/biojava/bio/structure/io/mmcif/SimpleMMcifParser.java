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
package org.biojava.bio.structure.io.mmcif;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.MMCIFFileReader;
import org.biojava.bio.structure.io.StructureIOFile;
import org.biojava.bio.structure.io.mmcif.model.AtomSite;
import org.biojava.bio.structure.io.mmcif.model.AuditAuthor;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava.bio.structure.io.mmcif.model.ChemCompAtom;
import org.biojava.bio.structure.io.mmcif.model.ChemCompBond;
import org.biojava.bio.structure.io.mmcif.model.ChemCompDescriptor;
import org.biojava.bio.structure.io.mmcif.model.DatabasePDBremark;
import org.biojava.bio.structure.io.mmcif.model.DatabasePDBrev;
import org.biojava.bio.structure.io.mmcif.model.Entity;
import org.biojava.bio.structure.io.mmcif.model.EntityPolySeq;
import org.biojava.bio.structure.io.mmcif.model.Exptl;
import org.biojava.bio.structure.io.mmcif.model.PdbxChemCompDescriptor;
import org.biojava.bio.structure.io.mmcif.model.PdbxChemCompIdentifier;
import org.biojava.bio.structure.io.mmcif.model.PdbxEntityNonPoly;
import org.biojava.bio.structure.io.mmcif.model.PdbxNonPolyScheme;
import org.biojava.bio.structure.io.mmcif.model.PdbxPolySeqScheme;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructOperList;
import org.biojava.bio.structure.io.mmcif.model.Refine;
import org.biojava.bio.structure.io.mmcif.model.Struct;
import org.biojava.bio.structure.io.mmcif.model.StructAsym;
import org.biojava.bio.structure.io.mmcif.model.StructKeywords;
import org.biojava.bio.structure.io.mmcif.model.StructRef;
import org.biojava.bio.structure.io.mmcif.model.StructRefSeq;
import org.biojava.bio.structure.jama.Matrix;


/** A simple mmCif file parser
 *
 * @author Andreas Prlic
 * @since 1.7
 * Usage:
 * <pre>
   		String file = "path/to/mmcif/file";
  		StructureIOFile pdbreader = new MMCIFFileReader();
		try {
			Structure s = pdbreader.getStructure(file);
			System.out.println(s);

			// you can convert it to a PDB file...
			System.out.println(s.toPDB());
		} catch (IOException e) {
			e.printStackTrace();
		}
 * </pre>
 * For more documentation see <a href="http://biojava.org/wiki/BioJava:CookBook#Protein_Structure">http://biojava.org/wiki/BioJava:CookBook#Protein_Structure</a>.
 */
public class SimpleMMcifParser implements MMcifParser {

	List<MMcifConsumer> consumers ;

	public static final String LOOP_END = "#";
	public static final String LOOP_START = "loop_";
	public static final String FIELD_LINE = "_";
	public static final String STRING_LIMIT = ";";

	private static final char s1 = '\'';
	private static final char s2 = '\"';

	Struct struct ;

	public static Logger logger =  Logger.getLogger("org.biojava.bio.structure");

	public SimpleMMcifParser(){
		consumers = new ArrayList<MMcifConsumer>();
		struct = null;
	}

	public void addMMcifConsumer(MMcifConsumer consumer) {
		consumers.add(consumer);

	}

	public void clearConsumers() {
		consumers.clear();

	}

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

	public void parse(InputStream inStream) throws IOException {
		parse(new BufferedReader(new InputStreamReader(inStream)));

	}

	public void parse(BufferedReader buf)
	throws IOException {

		triggerDocumentStart();


		// init container objects...
		struct = new Struct();
		String line = null;

		boolean inLoop = false;

		List<String> loopFields = new ArrayList<String>();
		List<String> lineData   = new ArrayList<String>();

		String category = null;


		// the first line is a data_PDBCODE line, test if this looks like a mmcif file
		line = buf.readLine();
		if (!line.startsWith("data_")){
			System.err.println("this does not look like a valid MMcif file! The first line should be data_1XYZ, but is " + line);
			triggerDocumentEnd();
			return;
		}

		while ( (line = buf.readLine ()) != null ){
			//System.out.println(inLoop + " " + line);


			if ( inLoop){

				if (line.startsWith(LOOP_END)){
					// reset all data
					inLoop = false;
					lineData.clear();
					category=null;
					loopFields.clear();
					continue;


				}

				if ( line.startsWith(FIELD_LINE)){
					// found another field.
					String txt = line.trim();
					//System.out.println("line: " + txt);
					if ( txt.indexOf('.') > -1){

						String[] spl = txt.split("\\.");
						//System.out.println(spl.length);
						category = spl[0];
						String attribute = spl[1];
						loopFields.add(attribute);
						if ( spl.length > 2){
							System.err.println("found nested attribute, not supported, yet!");
						}
					} else {
						category = txt;
					}


				} else {

					// in loop and we found a data line
					lineData = processLine(line, buf, loopFields.size());
					if ( lineData.size() != loopFields.size()){
						System.err.println("did not find enough data fields...");

					}

					endLineChecks(category, loopFields,lineData);

					lineData.clear();

				}

			} else {
				// not in loop

				if ( line.startsWith(LOOP_START)){
					loopFields.clear();
					inLoop = true;
					category=null;
					lineData.clear();
					continue;
				} else if (line.startsWith(LOOP_END)){
					inLoop = false;
					if ( category != null)
						endLineChecks(category, loopFields, lineData);
					category = null;
					loopFields.clear();
					lineData.clear();
				} else {
					// a boring normal line
					//System.out.println("boring data line: " + line + " " + inLoop + " " );
					List<String> data = processLine(line, buf, 2);
					//System.out.println("got a single line " + data);
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
						if (! line.startsWith("data_")){
							System.err.println("this does not look like a valid MMcif file! The first line should be data_1XYZ, but is " + line);
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


				}
			}
		}

		if (struct != null){
			triggerStructData(struct);
		}

		triggerDocumentEnd();

	}

	private List<String> processSingleLine(String line){
		//System.out.println("SS processSingleLine " + line);
		List<String> data = new ArrayList<String>();

		if ( line.trim().length() == 0){
			return data;
		}
		
		if ( line.trim().length() == 1){
			if ( line.startsWith(STRING_LIMIT))
				return data;
		}
		boolean inString = false;
		boolean inS1     = false;
		boolean inS2     = false;
		String word 	 = "";

		//System.out.println(line);
		for (int i=0; i< line.length(); i++ ){
			//System.out.println(word);
			Character c = line.charAt(i);

			Character nextC = null;
			if (i < line.length() - 1)
				nextC = line.charAt(i+1);
			if  (c == ' ') {

				if ( ! inString){
					if ( ! word.equals(""))
						data.add(word.trim());
					word = "";
				} else {
					// we are in a string, add the space
					word += c;
				}

			} else if (c == s1 )  {

				if ( inString){

					boolean wordEnd = false;
					if (! inS2) {
						if (nextC != null){
							//System.out.println("nextC: >"+nextC+"<");
							if ( Character.isWhitespace(nextC)){
								i++;
								wordEnd = true;
							}
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

				} else {
					// the beginning of a new string
					inString = true;
					inS1     = true;
				}
			} else if ( c == s2 ){
				if ( inString){

					boolean wordEnd = false;
					if (! inS1) {
						if (nextC != null){
							//System.out.println("nextC: >"+nextC+"<");
							if ( Character.isWhitespace(nextC)){
								i++;
								wordEnd = true;
							}
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
				} else {
					// the beginning of a new string
					inString = true;
					inS2     = true;
				}
			} else {
				word += c;
			}

		}
		if ( ! word.trim().equals(""))
			data.add(word);

		
		return data;

	}

	/** get the content of a cif entry
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

		String bigWord = null;

		while ( true ){

			if ( line.startsWith(STRING_LIMIT)){
				if (! inString){

					inString = true;
					if ( line.length() > 1)
						bigWord = line.substring(1);
					else
						bigWord = "";


				} else {
					// the end of a word
					lineData.add(bigWord);
					bigWord = null;
					inString = false;

				}
			} else {
				if ( inString )
					//TODO: make bigWord a stringbuffer...
					bigWord += (line);
				else {

					List<String> dat = processSingleLine(line);

					for (String d : dat){
						lineData.add(d);
					}
				}
			}

			//System.out.println("in process line : " + lineData.size() + " " + fieldLength);

			if ( lineData.size() > fieldLength){

				System.err.println("wrong data length ("+lineData.size()+
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


	private void endLineChecks(String category,List<String> loopFields, List<String> lineData ) throws IOException{


		/*System.out.println("parsed the following data: " +category + " fields: "+
				loopFields + " DATA: " +
				lineData);
		if (category.equals("_struct")){
			System.exit(0);
		}*/
		if ( loopFields.size() != lineData.size()){
			System.err.println("looks like we got a problem with nested string quote characters:");
			throw new IOException("data length ("+ lineData.size() +
					") != fields length ("+loopFields.size()+
					") category: " +category + " fields: "+
					loopFields + " DATA: " +
					lineData );
		}

		if ( category.equals("_entity")){

			Entity e =  (Entity) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.Entity",
					loopFields,lineData);
			triggerNewEntity(e);

		} else if ( category.equals("_struct")){

			struct =  (Struct) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.Struct",
					loopFields, lineData);

		} else if ( category.equals("_atom_site")){

			AtomSite a = (AtomSite) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.AtomSite",
					loopFields, lineData);
			triggerNewAtomSite(a);

		} else if ( category.equals("_database_PDB_rev")){
			DatabasePDBrev dbrev = (DatabasePDBrev) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.DatabasePDBrev",
					loopFields, lineData);

			triggerNewDatabasePDBrev(dbrev);

		} else if (  category.equals("_database_PDB_remark")){
			DatabasePDBremark remark = (DatabasePDBremark) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.DatabasePDBremark",
					loopFields, lineData);

			triggerNewDatabasePDBremark(remark);

		} else if ( category.equals("_exptl")){
			Exptl exptl  = (Exptl) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.Exptl",
					loopFields,lineData);

			triggerExptl(exptl);

		} else if ( category.equals("_struct_ref")){
			StructRef sref  = (StructRef) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.StructRef",
					loopFields,lineData);

			triggerNewStrucRef(sref);

		} else if ( category.equals("_struct_ref_seq")){
			StructRefSeq sref  = (StructRefSeq) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.StructRefSeq",
					loopFields,lineData);

			triggerNewStrucRefSeq(sref);
		} else if ( category.equals("_entity_poly_seq")){
			EntityPolySeq exptl  = (EntityPolySeq) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.EntityPolySeq",
					loopFields,lineData);

			triggerNewEntityPolySeq(exptl);

		} else if ( category.equals("_struct_asym")){
			StructAsym sasym  = (StructAsym) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.StructAsym",
					loopFields,lineData);

			triggerNewStructAsym(sasym);

		} else if ( category.equals("_pdbx_poly_seq_scheme")){
			PdbxPolySeqScheme ppss  = (PdbxPolySeqScheme) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.PdbxPolySeqScheme",
					loopFields,lineData);

			triggerNewPdbxPolySeqScheme(ppss);

		} else if ( category.equals("_pdbx_nonpoly_scheme")){
			PdbxNonPolyScheme ppss  = (PdbxNonPolyScheme) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.PdbxNonPolyScheme",
					loopFields,lineData);

			triggerNewPdbxNonPolyScheme(ppss);

		} else if ( category.equals("_pdbx_entity_nonpoly")){
			PdbxEntityNonPoly pen = (PdbxEntityNonPoly) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.PdbxEntityNonPoly",
					loopFields,lineData
			);
			triggerNewPdbxEntityNonPoly(pen);
		} else if ( category.equals("_struct_keywords")){
			StructKeywords kw = (StructKeywords)buildObject(
					"org.biojava.bio.structure.io.mmcif.model.StructKeywords",
					loopFields,lineData
			);
			triggerNewStructKeywords(kw);
		} else if (category.equals("_refine")){
			Refine r = (Refine)buildObject(
					"org.biojava.bio.structure.io.mmcif.model.Refine",
					loopFields,lineData
			);
			triggerNewRefine(r);
		} else if (category.equals("_chem_comp")){
			ChemComp c = (ChemComp)buildObject(
					"org.biojava.bio.structure.io.mmcif.model.ChemComp",
					loopFields, lineData
					);
			
			triggerNewChemComp(c);
		} else if (category.equals("_audit_author")) {
		   AuditAuthor aa = (AuditAuthor)buildObject(
		         "org.biojava.bio.structure.io.mmcif.model.AuditAuthor",
		         loopFields, lineData);
		      triggerNewAuditAuthor(aa);
		} else if (category.equals("_pdbx_chem_comp_descriptor")) {
			ChemCompDescriptor ccd = (ChemCompDescriptor) buildObject(
			         "org.biojava.bio.structure.io.mmcif.model.ChemCompDescriptor",
			         loopFields, lineData);
			triggerNewChemCompDescriptor(ccd);
		} else if (category.equals("_pdbx_struct_oper_list")) {
			/* PdbxStructOperList structOper = (PdbxStructOperList) buildObject(
					"org.biojava.bio.structure.io.mmcif.model.PdbxStructOperList",
			         loopFields, lineData); */
			
			// this guy is special since we convert to Matrices and shift vectors...
			PdbxStructOperList structOper = getPdbxStructOperList(loopFields,lineData);
			triggerNewPdbxStructOper(structOper);
			
			
		} else if (category.equals("_pdbx_struct_assembly")) {
			PdbxStructAssembly sa = (PdbxStructAssembly) buildObject(
			         "org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly",
			         loopFields, lineData);			
			triggerNewPdbxStructAssembly(sa);
			
		} else if (category.equals("_pdbx_struct_assembly_gen")) {
			PdbxStructAssemblyGen sa = (PdbxStructAssemblyGen) buildObject(
			         "org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGen",
			         loopFields, lineData);			
			triggerNewPdbxStructAssemblyGen(sa);
		} else if ( category.equals("_chem_comp_atom")){
			ChemCompAtom atom = (ChemCompAtom)buildObject(
					"org.biojava.bio.structure.io.mmcif.model.ChemCompAtom",
					loopFields,lineData);
			triggerNewChemCompAtom(atom);
			
		}else if ( category.equals("_chem_comp_bond")){
			ChemCompBond bond = (ChemCompBond)buildObject(
					"org.biojava.bio.structure.io.mmcif.model.ChemCompBond",
					loopFields,lineData);
			triggerNewChemCompBond(bond);
		} else if ( category.equals("_pdbx_chem_comp_identifier")){
			PdbxChemCompIdentifier id = (PdbxChemCompIdentifier)buildObject(
					"org.biojava.bio.structure.io.mmcif.model.PdbxChemCompIdentifier",
					loopFields,lineData);
			triggerNewPdbxChemCompIdentifier(id);
		} else if ( category.equals("_pdbx_chem_comp_descriptor")){
				PdbxChemCompDescriptor id = (PdbxChemCompDescriptor)buildObject(
						"org.biojava.bio.structure.io.mmcif.model.PdbxChemCompDescriptor",
						loopFields,lineData);
				triggerNewPdbxChemCompDescriptor(id);
						
		} else {
		
		

			// trigger a generic bean that can deal with all missing data types...
			triggerGeneric(category,loopFields,lineData);
		}


	}

	

//	@SuppressWarnings({ "rawtypes", "unchecked"})
//	private void setPair(Object o, List<String> lineData){
//		Class c = o.getClass();
//
//		if (lineData.size() == 2){
//			String key = lineData.get(0);
//			String val = lineData.get(1);
//
//			int dotPos = key.indexOf('.');
//
//			if ( dotPos > -1){
//				key = key.substring(dotPos+1,key.length());
//			}
//
//			String u = key.substring(0,1).toUpperCase();
//			try {
//				Method m = c.getMethod("set" + u + key.substring(1,key.length()) , String.class);
//				m.invoke(o,val);
//			}
//			catch (InvocationTargetException iex){
//				iex.printStackTrace();
//			}
//			catch (IllegalAccessException aex){
//				aex.printStackTrace();
//			}
//			catch( NoSuchMethodException nex){
//				if ( val.equals("?") || val.equals(".")) {
//					logger.info("trying to set field >" + key + "< in >"+ c.getName() + "<, but not found. Since value is >"+val+"<  most probably just ignore this.");
//				} else {
//					logger.warning("trying to set field >" + key + "< in >"+ c.getName() + "<, but not found! (value:" + val + ")");
//				}
//			}
//		} else {
//			System.err.println("trying to set key/value pair on object " +o.getClass().getName() + " but did not find in " + lineData);
//		}
//	}

	

	

	private PdbxStructOperList getPdbxStructOperList(List<String> loopFields,
			List<String> lineData) {
		PdbxStructOperList so = new PdbxStructOperList();
		
		//System.out.println(loopFields);
		//System.out.println(lineData);
		
		String id = lineData.get(loopFields.indexOf("id"));
		so.setId(id);
		so.setType(lineData.get(loopFields.indexOf("type")));
		Matrix matrix = new Matrix(3,3);
		for (int i = 1 ; i <=3 ; i++){
			for (int j =1 ; j <= 3 ; j++){
				String max = String.format("matrix[%d][%d]",j,i);
		
				String val = lineData.get(loopFields.indexOf(max));
				Double d = Double.parseDouble(val);
				 matrix.set(i-1,j-1,d);
			}
		}
		
		double[] coords =new double[3];
		
		for ( int i = 1; i <=3 ; i++){
			String v = String.format("vector[%d]",i);
			String val = lineData.get(loopFields.indexOf(v));
			Double d = Double.parseDouble(val);
			coords[i-1] = d;
		}
		
		so.setMatrix(matrix);
		so.setVector(coords);
		
		
		
		return so;
	}

	public void triggerNewPdbxStructOper(PdbxStructOperList structOper) {
		for(MMcifConsumer c : consumers){
			c.newPdbxStructOperList(structOper);
		}
		
	}

	@SuppressWarnings("rawtypes")
	private void setArray(Class c, Object o, String key, String val){


		// TODO: not implemented yet!
			//logger.info("Setting of array not implemented at the present for " + key + " " + val);
		/*
		int pos = key.indexOf("[");
		String varName = key.substring(0,pos);
		String u = varName.substring(0,1).toUpperCase();
		try {
			Method m = c.getMethod("set" + u + varName.substring(1,varName.length()) , String.class);
			m.invoke(o,val);
		} catch (Exception e){
			e.printStackTrace();
		}
		 */

	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	private Object buildObject(String className, List<String> loopFields, List<String> lineData) {
		Object o = null;
		try {
			// build up the Entity object from the line data...
			Class c = Class.forName(className);

			o = c.newInstance();


			int pos = -1 ;
			for (String key: loopFields){
				pos++;

				String val = lineData.get(pos);
				//System.out.println(key + " " + val);
				String u = key.substring(0,1).toUpperCase();

				try {
					Method m = c.getMethod("set" + u + key.substring(1,key.length()) , String.class);
					m.invoke(o,val);
				}
				catch( NoSuchMethodException nex){

					if (key.indexOf("[") > -1) {
						setArray(c,o,key,val);

					} else {
						System.err.println("Trying to set field " + key + " in "+ c.getName() +", but not found! (value:" + val + ")");
						//System.err.println(lineData);
					}
				}
			}
		} catch (InstantiationException eix){
			eix.printStackTrace();
		} catch (InvocationTargetException etx){
			etx.printStackTrace();
		} catch (IllegalAccessException eax){
			eax.printStackTrace();
		} catch (ClassNotFoundException ex){
			ex.printStackTrace();
		}
		return o;
	}

	public void triggerGeneric(String category, List<String> loopFields, List<String> lineData){
		for(MMcifConsumer c : consumers){
			c.newGenericData(category,loopFields, lineData);
		}
	}

	public void triggerNewEntity(Entity entity){
		for(MMcifConsumer c : consumers){
			c.newEntity(entity);
		}
	}

	public void triggerNewEntityPolySeq(EntityPolySeq epolseq){
		for(MMcifConsumer c : consumers){
			c.newEntityPolySeq(epolseq);
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
	private void triggerNewDatabasePDBrev(DatabasePDBrev dbrev){
		for(MMcifConsumer c : consumers){
			c.newDatabasePDBrev(dbrev);
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
}
