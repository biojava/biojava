/*
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
 * Created on 26.04.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.bio.structure.io;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.DBRef;
import org.biojava.bio.structure.Element;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.SSBond;
import org.biojava.bio.structure.Site;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava3.core.util.XMLWriter;


/** Methods to convert a structure object into different file formats.
 * @author Andreas Prlic
 * @since 1.4
 */
public class FileConvert {
	Structure structure ;

	boolean printConnections;

	// Locale should be english, e.g. in DE separator is "," -> PDB files have "." !
	static DecimalFormat d3 = (DecimalFormat)NumberFormat.getInstance(java.util.Locale.UK);
	static {
		d3.setMaximumIntegerDigits(3);
		d3.setMinimumFractionDigits(3);
		d3.setMaximumFractionDigits(3);
	}
	static DecimalFormat d2 = (DecimalFormat)NumberFormat.getInstance(java.util.Locale.UK);
	static {
		d2.setMaximumIntegerDigits(2);
		d2.setMinimumFractionDigits(2);
		d2.setMaximumFractionDigits(2);
	}

	private static final String newline = System.getProperty("line.separator");

	/**
	 * Constructs a FileConvert object.
	 *
	 * @param struc  a Structure object
	 */
	public FileConvert(Structure struc) {
		structure = struc ;
		printConnections = true;
	}

	/** align a string to the right
	 * length is the total length the new string should take, inlcuding spaces on the left
	 * incredible that this tool is missing in java !!!
	 */
	private static String alignRight(String input, int length){


		int n = input.length();
		if ( n >= length)
			return input;

		String spaces = "                           " ;
		int diff = length - n ;
		StringBuffer s = new StringBuffer();

		s.append(spaces.substring(0,diff));
		s.append(input);

		return s.toString();
	}

	private static String alignLeft(String input, int length){
		if (input.length() >= length) {
			return input;
		}

		String spaces = "                           " ;
		input += spaces.substring(0, length - input.length() );
		return input;

	}

	/** returns if the Connections should be added
	 * default is true;
	 * @return if the printConnections flag is set
	 */
	public boolean doPrintConnections() {
		return printConnections;
	}

	/** enable/disable printing of connections
	 * connections are sometimes buggy in PDB files
	 * so there are some cases where one might turn this off.
	 * @param printConnections
	 */
	public void setPrintConnections(boolean printConnections) {
		this.printConnections = printConnections;
	}

	/** prints the connections in PDB style
	 *
	 * Thanks to Tamas Horvath for this one
	 */
	private String printPDBConnections(){


		StringBuffer str = new StringBuffer();

		List<Map<String, Integer>> cons = structure.getConnections();
		for (int cnr = 0; cnr<cons.size();cnr++){
			Map<String,Integer> con =  cons.get(cnr);
			Integer as = (Integer) con.get("atomserial");

			String atomserial =  "";

			String bond1 = "";
			String bond2 = "";
			String bond3 = "";
			String bond4 = "";
			String hyd1  = "";
			String hyd2  = "";
			String salt1 = "";
			String hyd3  = "";
			String hyd4  = "";
			String salt2 = "";



			if (con.containsKey("bond1"))  bond1 = con.get("bond1").toString();
			if (con.containsKey("bond2"))  bond2 = con.get("bond2").toString();
			if (con.containsKey("bond3"))  bond3 = con.get("bond3").toString();
			if (con.containsKey("bond4"))  bond4 = con.get("bond4").toString();
			if (con.containsKey("hyd1"))   hyd1  = con.get("hyd1").toString();
			if (con.containsKey("hyd2"))   hyd2  = con.get("hyd2").toString();
			if (con.containsKey("salt1"))  salt1 = con.get("salt1").toString();
			if (con.containsKey("hyd3"))   hyd3  = con.get("hyd3").toString();
			if (con.containsKey("hyd4"))   hyd4  = con.get("hyd4").toString();
			if (con.containsKey("salt2"))  salt2 = con.get("salt2").toString();

			atomserial = alignRight(""+as,5) ;
			bond1      = alignRight(bond1,5) ;
			bond2      = alignRight(bond2,5) ;
			bond3      = alignRight(bond3,5) ;
			bond4      = alignRight(bond4,5) ;
			hyd1       = alignRight(hyd1,5)  ;
			hyd2       = alignRight(hyd2,5)  ;
			salt1      = alignRight(salt1,5) ;
			hyd3       = alignRight(hyd3,5)  ;
			hyd4       = alignRight(hyd4,5)  ;
			salt2      = alignRight(salt2,5) ;

			String connectLine = "CONECT" + atomserial + bond1 + bond2 + bond3 +
			bond4 + hyd1 + hyd2 + salt1 + hyd3 + hyd4 + salt2;

            str.append(connectLine).append(newline);
		}
		return str.toString();
	}

	/** Convert a structure into a PDB file.
	 * @return a String representing a PDB file.
	 */
	public String toPDB() {


		StringBuffer str = new StringBuffer();
		//int i = 0 ;



		// TODO: print all the PDB header informaton in PDB style
		// some objects (PDBHeader, Compound) are still missing
		//

		PDBHeader header = structure.getPDBHeader();
		header.toPDB(str);


		//REMARK 800
		if (!structure.getSites().isEmpty()) {
            str.append("REMARK 800                                                                      ").append(newline);
            str.append("REMARK 800 SITE                                                                 ").append(newline);
			for (Site site : structure.getSites()) {
				site.remark800toPDB(str);
			}
		}
		//DBREF
		for (DBRef dbref : structure.getDBRefs()){
			dbref.toPDB(str);
			str.append(newline);
		}
		//SSBOND
		for (SSBond ssbond : structure.getSSBonds()){
			ssbond.toPDB(str);
			str.append(newline);
		}
		//SITE
		for (Site site : structure.getSites()) {
			try {
				site.toPDB(str);
			} catch (Exception e){
				e.printStackTrace();
			}
		}

		//
		// print the atom records
		//

		// do for all models
		int nrModels = structure.nrModels() ;
		if ( structure.isNmr()) {
			str.append("EXPDTA    NMR, "+ nrModels+" STRUCTURES"+newline) ;
		}
		for (int m = 0 ; m < nrModels ; m++) {
			List<Chain> model = structure.getModel(m);
			// todo support NMR structures ...
			if ( structure.isNmr()) {
				str.append("MODEL      " + (m+1)+ newline);
			}
			// do for all chains
			int nrChains = model.size();
			for ( int c =0; c<nrChains;c++) {
				Chain  chain   = model.get(c);
				//String chainID = chain.getChainID();
				//if ( chainID.equals(DEFAULTCHAIN) ) chainID = " ";
				// do for all groups
				int nrGroups = chain.getAtomLength();
				for ( int h=0; h<nrGroups;h++){

					Group g= chain.getAtomGroup(h);

					
					toPDB(g,str);
					
					
				}
			}

			if ( structure.isNmr()) {
                str.append("ENDMDL").append(newline);
			}



		}

		if ( doPrintConnections() )
			str.append(printPDBConnections());

		return str.toString() ;
	}

	private static void toPDB(Group g, StringBuffer str) {
		// iterate over all atoms ...
		// format output ...
		int groupsize  = g.size();

		for ( int atompos = 0 ; atompos < groupsize; atompos++) {
			Atom a = null ;
			try {
				a = g.getAtom(atompos);
			} catch ( StructureException e) {
				System.err.println(e);
				continue ;
			}

			toPDB(a, str);


			//line = record + serial + " " + fullname +altLoc
			//+ leftResName + " " + chainID + resseq
			//+ "   " + x+y+z
			//+ occupancy + tempfactor;
			//str.append(line + newline);
			//System.out.println(line);
		}
		if ( g.hasAltLoc()){
			for (Group alt : g.getAltLocs() ) {
				toPDB(alt,str);
			}
		}
		
	}

	/** Prints the content of an Atom object as a PDB formatted line.
	 * 
	 * @param a
	 * @return
	 */
	public static String toPDB(Atom a){
		StringBuffer w = new StringBuffer();

		toPDB(a,w);

		return w.toString();

	}
	
	
	/** Convert a CHain object to PDB representation
	 * 
	 * @param chain
	 * @return
	 */
	public static String toPDB(Chain chain){
		StringBuffer w = new StringBuffer();
		int nrGroups = chain.getAtomLength();
		
		for ( int h=0; h<nrGroups;h++){

			Group g= chain.getAtomGroup(h);

			
			toPDB(g,w);
			
			
		}
		
		return w.toString();
	}
	
	/** Convert a Group object to PDB representation
	 * 
	 * @param g
	 * @return
	 */
	public static String toPDB(Group g){
		StringBuffer w = new StringBuffer();
		toPDB(g,w);
		return w.toString();
	}

	/**
	 * print ATOM record in the following syntax
	<pre>
    ATOM      1  N   ASP A  15     110.964  24.941  59.191  1.00 83.44           N
*
COLUMNS        DATA TYPE       FIELD         DEFINITION
---------------------------------------------------------------------------------
1 -  6        Record name     "ATOM  "
7 - 11        Integer         serial        Atom serial number.
13 - 16        Atom            name          Atom name.
17             Character       altLoc        Alternate location indicator.
18 - 20        Residue name    resName       Residue name.
22             Character       chainID       Chain identifier.
23 - 26        Integer         resSeq        Residue sequence number.
27             AChar           iCode         Code for insertion of residues.
31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
Angstroms.
39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
Angstroms.
47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
Angstroms.
55 - 60        Real(6.2)       occupancy     Occupancy.
61 - 66        Real(6.2)       tempFactor    Temperature factor.
73 - 76        LString(4)      segID         Segment identifier, left-justified.
77 - 78        LString(2)      element       Element symbol, right-justified.
79 - 80        LString(2)      charge        Charge on the atom.
</pre>
*/
	public static void toPDB(Atom a, StringBuffer str) {

		Group g = a.getGroup();
		Chain c = g.getChain();
		String chainID = c.getChainID();

		String type = g.getType() ;

		String record = "" ;
		if ( type.equals("hetatm") ) {
			record = "HETATM";
		} else {
			record = "ATOM  ";
		}


		// format output ...
		//int groupsize  = g.size();
		String resName = g.getPDBName(); 
		String pdbcode = g.getResidueNumber().toString();
		//String line    = "" ;


		int    seri       = a.getPDBserial()        ;
		String serial     = alignRight(""+seri,5)   ;
		String fullname   = a.getFullName()         ;

		// System.out.println(" fullname: " + fullname + " : " + a.getAltLoc() + " : " + pdbcode);

		Character  altLoc = a.getAltLoc()           ;
		String resseq = "" ;
		if ( hasInsertionCode(pdbcode) )
			resseq     = alignRight(""+pdbcode,5);
		else
			resseq     = alignRight(""+pdbcode,4)+" ";
		String x          = alignRight(""+d3.format(a.getX()),8);
		String y          = alignRight(""+d3.format(a.getY()),8);
		String z          = alignRight(""+d3.format(a.getZ()),8);
		String occupancy  = alignRight(""+d2.format(a.getOccupancy()),6) ;
		String tempfactor = alignRight(""+d2.format(a.getTempFactor()),6);

		//System.out.println("fullname,zise:" + fullname + " " + fullname.length());

		String leftResName = alignLeft(resName,3);

		StringBuffer s = new StringBuffer();
		s.append(record);
		s.append(serial);
		s.append(" ");
		s.append(fullname);
		s.append(altLoc);
		s.append(leftResName);
		s.append(" ");
		s.append(chainID);
		s.append(resseq);
		s.append("   ");
		s.append(x);
		s.append(y);
		s.append(z);
		s.append(occupancy);
		s.append(tempfactor);

		Element e = a.getElement();
		
		String eString = e.toString().toUpperCase();
		
		if ( e.equals(Element.R)) {
			eString = "X";
		}
		str.append(String.format("%-76s%2s", s.toString(),eString));
		str.append(newline);

	}

	/** test if pdbserial has an insertion code */
	private static boolean hasInsertionCode(String pdbserial) {
		try {
			Integer.parseInt(pdbserial) ;
		} catch (NumberFormatException e) {
			return true ;
		}
		return false ;
	}


	/** convert a protein Structure to a DAS Structure XML response .
	 * @param xw  a XMLWriter object
	 * @throws IOException ...
	 *
	 */
	@SuppressWarnings("deprecation")
	public void toDASStructure(XMLWriter xw)
	throws IOException
	{

		/*xmlns="http://www.sanger.ac.uk/xml/das/2004/06/17/dasalignment.xsd" xmlns:align="http://www.sanger.ac.uk/xml/das/2004/06/17/alignment.xsd" xmlns:xsd="http://www.w3.org/2001/XMLSchema-instance" xsd:schemaLocation="http://www.sanger.ac.uk/xml/das/2004/06/17/dasalignment.xsd http://www.sanger.ac.uk/xml/das//2004/06/17/dasalignment.xsd"*/

		if ( structure == null){
			System.err.println("can not convert structure null");
			return;
		}

		Map<String,Object> header = structure.getHeader();

		xw.openTag("object");
		xw.attribute("dbAccessionId",structure.getPDBCode());
		xw.attribute("intObjectId"  ,structure.getPDBCode());
		// missing modification date
		String modificationDate = (String)header.get("modDate") ;
		xw.attribute("objectVersion",modificationDate);
		xw.attribute("type","protein structure");
		xw.attribute("dbSource","PDB");
		xw.attribute("dbVersion","20070116");
		xw.attribute("dbCoordSys","PDBresnum,Protein Structure");

		// do we need object details ???
		xw.closeTag("object");


		// do for all models
		for (int modelnr = 0;modelnr<structure.nrModels();modelnr++){

			// do for all chains:
			for (int chainnr = 0;chainnr<structure.size(modelnr);chainnr++){
				Chain chain = (Chain)structure.getChain(modelnr,chainnr);
				xw.openTag("chain");
				xw.attribute("id",chain.getChainID());
				xw.attribute("SwissprotId",chain.getSwissprotId() );
				if (structure.isNmr()){
					xw.attribute("model",Integer.toString(modelnr+1));
				}

				//do for all groups:
				for (int groupnr =0;
				groupnr<chain.getAtomLength()
				;groupnr++){
					Group gr = chain.getAtomGroup(groupnr);
					xw.openTag("group");
					xw.attribute("name",gr.getPDBName());
					xw.attribute("type",gr.getType());
					xw.attribute("groupID",gr.getResidueNumber().toString());


					// do for all atoms:
					//Atom[] atoms  = gr.getAtoms();
					List<Atom> atoms =  gr.getAtoms();
					for (int atomnr=0;atomnr<atoms.size();atomnr++){
						Atom atom = (Atom)atoms.get(atomnr);
						xw.openTag("atom");
						xw.attribute("atomID",Integer.toString(atom.getPDBserial()));
						xw.attribute("atomName",atom.getFullName());
						xw.attribute("x",Double.toString(atom.getX()));
						xw.attribute("y",Double.toString(atom.getY()));
						xw.attribute("z",Double.toString(atom.getZ()));
						xw.closeTag("atom");
					}
					xw.closeTag("group") ;
				}

				xw.closeTag("chain");
			}
		}


		if ( doPrintConnections() ) {
			// do connectivity for all chains:

			List<Map<String,Integer>> cons = structure.getConnections();
			for (int cnr = 0; cnr<cons.size();cnr++){


				/*
                 the HashMap for a single CONECT line contains the following fields:
                 <ul>
                 <li>atomserial (mandatory) : Atom serial number
                 <li>bond1 .. bond4 (optional): Serial number of bonded atom
                 <li>hydrogen1 .. hydrogen4 (optional):Serial number of hydrogen bonded atom
                 <li>salt1 .. salt2 (optional): Serial number of salt bridged atom
                 </ul>
				 */

				Map<String, Integer> con = (Map<String, Integer>)cons.get(cnr);
				Integer as = (Integer)con.get("atomserial");
				int atomserial = as.intValue();


				List<Integer> atomids = new ArrayList<Integer>() ;

				// test salt and hydrogen first //
				if (con.containsKey("salt1")) atomids.add(con.get("salt1"));
				if (con.containsKey("salt2")) atomids.add(con.get("salt2"));

				if (atomids.size()!=0){
					addConnection(xw,"salt",atomserial,atomids);
					atomids = new ArrayList<Integer>() ;
				}
				if (con.containsKey("hydrogen1")) atomids.add(con.get("hydrogen1"));
				if (con.containsKey("hydrogen2")) atomids.add(con.get("hydrogen2"));
				if (con.containsKey("hydrogen3")) atomids.add(con.get("hydrogen3"));
				if (con.containsKey("hydrogen4")) atomids.add(con.get("hydrogen4"));
				if (atomids.size()!=0){
					addConnection(xw,"hydrogen",atomserial,atomids);
					atomids = new ArrayList<Integer>() ;
				}

				if (con.containsKey("bond1")) atomids.add(con.get("bond1"));
				if (con.containsKey("bond2")) atomids.add(con.get("bond2"));
				if (con.containsKey("bond3")) atomids.add(con.get("bond3"));
				if (con.containsKey("bond4")) atomids.add(con.get("bond4"));

				if (atomids.size()!=0){
					addConnection(xw,"bond",atomserial,atomids);
				}
			}
		}
	}

	private void addConnection(XMLWriter xw,String connType, int atomserial, List<Integer> atomids){
		try{
			xw.openTag("connect");
			xw.attribute("atomSerial",Integer.toString(atomserial));
			xw.attribute("type",connType);
			for (int i=0;i<atomids.size();i++){
				Integer atomid = atomids.get(i);
				if ( atomid == null)
					continue;
				int aid = atomid.intValue();
				xw.openTag("atomID");
				xw.attribute("atomID",Integer.toString(aid));
				xw.closeTag("atomID");
			}
			xw.closeTag("connect");
		} catch( Exception e) {
			e.printStackTrace();
		}
	}


}
